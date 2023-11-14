import h5py
import numpy as np
import os
import time
import compartment
from common import \
    gk, gna, gcl, ghco3, gh, p_atpase, p_kcc2, h2co3_i, \
    pw, vw, RTF, cm, F, kf


class Simulator:

    def __init__(self, sim_type="New", new_file_name="", old_file_name=""):
        """
        Initializing a new instance of the simulator
        @param file_name: string of the file name
        """

        self.sim_type = sim_type
        self.new_file_name = new_file_name
        self.old_file_name = old_file_name

        # Creating a new HDF file
        current_folder = os.getcwd()
        parent_folder = os.path.dirname(os.path.abspath(current_folder))
        parent_folder2 = os.path.dirname(os.path.abspath(parent_folder)) #Note that HDF5 files are saved in a folder directly outside the repo as these files can be large.
        hdf5_folder = os.path.join(parent_folder2, "HDF5_files")
        # Create the "HDF5_files" folder if it doesn't exist
        if not os.path.exists(hdf5_folder):
            os.makedirs(hdf5_folder)
        self.new_file_name = os.path.join(hdf5_folder, self.new_file_name)
        self.old_file_name = os.path.join(hdf5_folder, self.old_file_name)

        self.intra = None  # Intracellular compartment object
        self.extra = None  # Extracellular compartment object

        self.sim_start_t = None
        self.dt = None
        self.total_t = None
        self.total_intervals = None
        self.total_steps = None
        self.interval_num = 0
        self.run_t = 0
        self.step_num = 0

        self.FinvC = F / cm

        self.infinite_bath = True
        self.KCC2_on = True
        self.dynamic_ATPase = True

        self.p_ATPase = p_atpase  # ATPase pump rate
        self.j_ATPase = p_atpase * (14e-3 / 145e-3) ** 3  # ATPase flux default
        self.p_kcc2 = p_kcc2  # KCC2 pump rate
        self.j_kcc2 = None

        self.z_change_on = False
        self.z_change_params = None
        self.KCC2_change_on = False
        self.KCC2_params = None

        self.syn_on = False
        self.syn_params = None
        self.cl_syn,self.hco3_syn = 0,0

        #self.kf = 10 ** (6)  ## Reduce?
        self.kf = kf
        self.kr = kf*5
        #self.kr = (-5e-11 + self.kf * 0.0076) / (0.01 * 63e-5)
        #self.kr = (-5e-11 + self.kf * 0.0016) / (0.01 * 63e-9)
        print("k_r:" +str(self.kr))
        # self.kr = 10 ** (12.4)
        self.k_na_h = 5*self.j_ATPase
        self.h_imbalance = 0


        self.g_extra = 0
        self.current_na = 0

        if sim_type == "Extend":
            self.load_old_HDFfile()

    def load_old_HDFfile(self):

        with h5py.File(self.old_file_name, mode='r') as hdf_old:
            C = hdf_old.get('COMPARTMENTS')
            C_group_arr = []
            comp_names_arr = list(C.keys())
            T = hdf_old.get('TIMING')
            total_t = T.get('TOTAL_T')[()]
            intervals = T.get('INTERVALS')[()]
            dt = T.get("DT")[()]
            total_steps = round(total_t / dt)
            interval_step = total_steps / intervals
            interval_arr = [int(interval_step * i) for i in range(intervals)]

            comp_arr_lastvalues = []  # saves all the last values of the base experiment

            # Looping through old compartments and saving last dataset values

            for _ in comp_names_arr:
                comp = C.get(_)
                steplist = list(comp.keys())
                rev_interval_arr = interval_arr[::-1]
                for j in rev_interval_arr:
                    interval = j
                    if str(interval) in str(steplist):
                        last_time_point = str(interval)
                        break
                last_dataset = comp.get(last_time_point)
                comp_arr_lastvalues.append(list(last_dataset))

        hdf_old.close()

        for i in range(len(comp_names_arr)):
            last_dataset = comp_arr_lastvalues[i]
            last_dataset[0] = 0
            ## creating a new compartment object and adding to the simulator
            rad = last_dataset[1]
            length = 25e-5  # the length is not saved in the HDF5 file

            comp = compartment.Compartment(compartment_name=comp_names_arr[i], radius=rad, length=length)
            comp.set_ion_properties(na_i=last_dataset[3], k_i=last_dataset[4], cl_i=last_dataset[5],
                                    hco3_i=last_dataset[6], h_i=last_dataset[7],
                                    x_i=last_dataset[8], z_i=last_dataset[9])
            comp.v = last_dataset[10]

            if "intracellular" in comp.name.lower():
                self.intra = comp
            elif "extracellular" in comp.name.lower():
                self.extra = comp

    def set_timing(self, total_t=0, dt=0, total_intervals=0):
        """
        @param total_t: float, total simulation time (in seconds)
        @param dt: float, time step of the simulation(in seconds) this determines how often calculations are repeated
             default is set to 1e-6 seconds.
        @param intervals: integer, determines how often to write simulation results to the HDF file
        :param self:
        """
        self.total_t = total_t
        self.dt = dt
        self.total_intervals = total_intervals

        # total number of iterations that will occur in the simulation
        self.total_steps = self.total_t / self.dt

        # percentage of simulation run time to display(output)
        self.status_percentages = (0.001, 0.005, 0.01, 0.1, 0.25, 0.5, 0.75, 1)
        # array of the iteration numbers for each output interval
        self.status_steps = [round(self.status_percentages[a] * self.total_steps, 0) for a in
                             range(len(self.status_percentages))]
        self.status_steps = tuple(self.status_steps)

        # interval step determine at which iteration number the sim will be saved
        self.interval_step = self.total_steps / self.total_intervals
        # interval array contains the iteration numbers where sim will be saved
        self.interval_arr = [round(self.interval_step * (i)) for i in range(self.total_intervals + 1)]
        self.interval_arr = tuple(self.interval_arr)

        self.intra.dt = self.dt
        self.extra.dt = self.dt

    def set_kcc2_off(self):
        self.p_kcc2 = 0


    def set_intracellular_properties(self, compartment_params={}):
        self.intra = compartment.Compartment(compartment_name=compartment_params["name"],
                                             radius=compartment_params["radius"],
                                             length=compartment_params["length"])

        self.intra.set_ion_properties(na_i=compartment_params["na"],
                                      k_i=compartment_params["k"],
                                      cl_i=compartment_params["cl"],
                                      hco3_i=compartment_params["hco3"],
                                      h_i=compartment_params["h"],
                                      x_i=compartment_params["X"],
                                      z_i=compartment_params["z"])

    def set_extracellular_properties(self, compartment_params={}):
        self.extra = compartment.Compartment(compartment_name=compartment_params["name"],
                                             radius=compartment_params["radius"],
                                             length=compartment_params["length"])

        self.extra.set_ion_properties(na_i=compartment_params["na"],
                                      k_i=compartment_params["k"],
                                      cl_i=compartment_params["cl"],
                                      hco3_i=compartment_params["hco3"],
                                      h_i=compartment_params["h"],
                                      x_i=compartment_params["X"],
                                      z_i=compartment_params["z"])

    def set_z_change(self, start_t=0, end_t=0,
                     z_change_amount=-0.2, adjust_cl=True, target_vm=None):
        self.z_change_on = True
        self.z_change_params = {"start_t": start_t,
                                "end_t": end_t,
                                "z_change_amount": z_change_amount,
                                "z_change_delta": (z_change_amount * self.dt) / (end_t - start_t),
                                "adjust_cl": adjust_cl,
                                "target_vm": target_vm}

    def set_KCC2_change(self, start_t=0, end_t=10, final_KCC2_value=0):
        self.KCC2_change_on = True
        KCC2_change = final_KCC2_value - self.p_kcc2
        self.KCC2_change_params = {"start_t": start_t,
                                   "end_t": end_t,
                                   "final_KCC2_value": final_KCC2_value,
                                   "KCC2_change_delta": (KCC2_change * self.dt) / (end_t - start_t)}

    def add_synapse(self, start_t=0, tau=250 * 1e-3, max_g=1e-9):
        """
            function to add a synapse to a particular compartment
            @param start_t: float,start time for synaptic input
            @param duration: float,duration of synaptic input
            @param max_neurotransmitter: float,max neurotransmitter concentration in moles/liter
            @param synapse_conductance: float, conductance of the synapse channel in Siemens, default is 1nS
            """
        self.syn_on = True
        self.syn_params = {"start_t": start_t,
                           "tau": tau,
                           "max_g": max_g}

    def synapse_step(self):
        # Model of GABAergic input - "Alpha function"
        # Moved away from a receptor binding model towards a simpler alpha function of modelling gaba
        # i_cl = g(Vm-ECl) * (gaba_fraction)
        # i_hco3 = g(Vm-EHCO3) * (1 - gaba_fraction)
        # g = conductance of the synapse at a particular time point
        # g = gmax*(t-onset) / (tau * e[-(t-onset-tau)/tau]
        # gmax = maximum conductance
        # tau = time to peak of the conductance change (e.g 0.1ms)

        syn_g = self.syn_params["max_g"] * (
                (self.run_t - self.syn_params["start_t"]) / self.syn_params["tau"]) * np.exp(
            -(self.run_t - self.syn_params["start_t"] - self.syn_params["tau"]) / self.syn_params["tau"])

        # GABA Fraction
        gaba_fraction = (self.intra.E_hco3 - self.intra.E_gaba) / (
                self.intra.E_hco3 - self.intra.E_cl)  # ref. chris paper

        # Chloride component
        I_cl = syn_g * (self.intra.v - self.intra.E_cl) * gaba_fraction
        I_cl = I_cl / F  # converting coulomb to mol
        I_cl = I_cl * self.dt  # getting the mol input for the timestep
        cl_change = I_cl / self.intra.w  # calculating concentration change
        self.cl_syn = cl_change

        # Bicarb component
        I_hco3 = syn_g * (self.intra.v - self.intra.E_hco3) * (1 - gaba_fraction)
        I_hco3 = I_hco3 / F  # converting coulomb to mol
        I_hco3 = I_hco3 * self.dt  # getting the mol input for the timestep
        hco3_change = I_hco3 / self.intra.w  # calculating concentration change
        self.hco3_syn = hco3_change


    def calc_voltages(self):
        intracellular_charge = self.intra.na_i + self.intra.k_i + self.intra.h_i + \
                               (self.intra.z_i * self.intra.x_i) - (self.intra.cl_i + self.intra.hco3_i)

        extracellular_charge = self.extra.na_i + self.extra.k_i + self.intra.h_i + \
                               (self.extra.z_i * self.extra.x_i) - (self.extra.cl_i + self.extra.hco3_i)

        self.intra.v = self.FinvC * (
                (intracellular_charge * self.intra.w) - (extracellular_charge * self.extra.w)) / self.intra.sa

        # Nernst equations
        self.intra.E_k = -1 * RTF * np.log(self.intra.k_i / self.extra.k_i)
        self.intra.E_na = -1 * RTF * np.log(self.intra.na_i / self.extra.na_i)
        self.intra.E_cl = RTF * np.log(self.intra.cl_i / self.extra.cl_i)
        self.intra.E_hco3 = RTF * np.log(self.intra.hco3_i / self.extra.hco3_i)
        self.intra.E_h = -1 * RTF * np.log(self.intra.h_i / self.extra.h_i)
        
        # GHK equation used to calculate GABA reversal potential
        numerator = 4 / 5 * self.intra.cl_i + 1 / 5 * self.intra.hco3_i
        denominator = 4 / 5 * self.extra.cl_i + 1 / 5 * self.extra.hco3_i
        self.intra.E_gaba = RTF * np.log(numerator / denominator)

    def sim_step(self):
        # Part 1: Calculate voltages
        self.calc_voltages()

        # Part 2: Calculate pump rates / channel or synapse conductances

        # 2.1: ATPase Rate
        if self.dynamic_ATPase:
            self.j_ATPase = self.p_ATPase * (self.intra.na_i / self.extra.na_i) ** 3
        # 2.2: KCC2 Rate
        if self.KCC2_change_on:
            if self.KCC2_change_params["start_t"] < self.run_t < self.KCC2_change_params["end_t"]:
                self.p_kcc2 += self.KCC2_change_params["KCC2_change_delta"]
            if self.run_t > self.KCC2_change_params["end_t"]:
                self.p_kcc2 = self.KCC2_change_params["final_KCC2_value"]
        self.j_kcc2 = self.p_kcc2 * (self.intra.E_k - self.intra.E_cl)

        # 2.3: Z change
        if self.z_change_on:
            if self.z_change_params["start_t"] < self.run_t < self.z_change_params["end_t"]:
                if (self.z_change_params["target_vm"] is None) or (self.intra.v < self.z_change_params["target_vm"]):
                    self.intra.z_i += self.z_change_params["z_change_delta"]
                    if self.z_change_params["adjust_cl"]:
                        self.intra.cl_i = self.intra.na_i + self.intra.k_i + self.intra.h_i + (self.intra.x_i * self.intra.z_i)

        # 2.4: Synapse step
        if self.syn_on:
            if self.run_t > self.syn_params["start_t"]:
                self.synapse_step()

        # 2.5: External current

        # Part 3: Calculate concentration changes

        # 3.1: Intracellular ions
        #Cl-
        d_cl_leak = + self.dt * self.intra.sa / self.intra.w * (gcl + self.g_extra) * (
                self.intra.v + RTF * np.log(self.extra.cl_i / self.intra.cl_i))
        d_cl_kcc2 = + self.dt * self.intra.sa / self.intra.w * self.j_kcc2
        self.intra.d_cl_i = d_cl_leak + d_cl_kcc2 + self.cl_syn

        #HCO3-
        d_hco3_leak = + self.dt * self.intra.sa / self.intra.w * (ghco3 + self.g_extra) * (
                self.intra.v + RTF * np.log(self.extra.hco3_i / self.intra.hco3_i))
        d_hco3_forwardRx = self.dt * self.kf * h2co3_i
        d_hco3_reverseRx = self.dt * self.kr * self.intra.hco3_i * self.intra.h_i

        self.intra.d_hco3_i = d_hco3_leak + (d_hco3_forwardRx - d_hco3_reverseRx) + self.hco3_syn

        # should convert to mols of HCO3 --> mols of H+ --> mols of Na+ : this is the same

        #Na+
        d_na_leak = - (self.dt * self.intra.sa / self.intra.w) * (gna + self.g_extra) * (
                self.intra.v + RTF * np.log(self.intra.na_i / self.extra.na_i))
        d_na_atpase = - self.dt * self.intra.sa / self.intra.w * (+3 * self.j_ATPase)
        d_na_current = self.current_na # only if external current is added
        #self.h_imbalance += d_hco3_forwardRx - d_hco3_reverseRx
        #self.h_imbalance = d_hco3_forwardRx - d_hco3_reverseRx
        d_na_NA_H_exchange = self.k_na_h*(self.intra.E_na - self.intra.E_h)
        #d_na_NA_H_exchange = self.k_na_h * self.h_imbalance # The net H+ produced by the reaction is replaced by Na+
        #self.h_imbalance -= d_na_NA_H_exchange
        self.intra.d_na_i = d_na_leak + d_na_atpase + d_na_current + d_na_NA_H_exchange

        #H+
        d_h_leak = - (self.dt * self.intra.sa / self.intra.w) * (gh + self.g_extra) * (
                self.intra.v + RTF * np.log(self.intra.h_i / self.extra.h_i))
        self.intra.d_h_i = d_h_leak + (d_hco3_forwardRx - d_hco3_reverseRx) - d_na_NA_H_exchange

        #K+
        d_k_leak = - self.dt * self.intra.sa / self.intra.w * (gk + self.g_extra) * (
                self.intra.v + RTF * np.log(self.intra.k_i / self.extra.k_i))
        d_k_atpase = - self.dt * self.intra.sa / self.intra.w * (- 2 * self.j_ATPase)
        d_k_kcc2 = - self.dt * self.intra.sa / self.intra.w * (- self.j_kcc2)
        self.intra.d_k_i = d_k_leak + d_k_atpase + d_k_kcc2


        self.intra.na_i = self.intra.na_i + self.intra.d_na_i
        self.intra.k_i = self.intra.k_i + self.intra.d_k_i
        self.intra.cl_i = self.intra.cl_i + self.intra.d_cl_i
        self.intra.hco3_i = self.intra.hco3_i + self.intra.d_hco3_i
        self.intra.h_i = self.intra.h_i + self.intra.d_h_i

        # 3.2: Extracellular ions (update this)
        if not self.infinite_bath:
            self.extra.na_i = self.extra.na_i - self.intra.d_na_i * self.intra.w / self.extra.w
            self.extra.k_i = self.extra.k_i - self.intra.d_k_i * self.intra.w / self.extra.w
            self.extra.cl_i = self.extra.cl_i - self.intra.d_cl_i * self.intra.w / self.extra.w

        # Part 4: Update cell volumes

        # 4.1: Calculate intracellular volume change
        intra_osm = self.intra.na_i + self.intra.k_i + self.intra.cl_i + self.intra.x_i + self.intra.hco3_i + self.intra.h_i
        extra_osm = self.extra.na_i + self.extra.k_i + self.extra.cl_i + self.extra.x_i + self.extra.hco3_i + self.intra.h_i
        dw = self.dt * (vw * pw * self.intra.sa * (intra_osm - extra_osm))
        intra_w2 = self.intra.w + dw

        # 4.2: Update intracellular ion concentrations
        self.intra.na_i = self.intra.na_i * self.intra.w / intra_w2
        self.intra.k_i = self.intra.k_i * self.intra.w / intra_w2
        self.intra.cl_i = self.intra.cl_i * self.intra.w / intra_w2
        self.intra.x_i = self.intra.x_i * self.intra.w / intra_w2
        self.intra.hco3_i = self.intra.hco3_i * self.intra.w / intra_w2
        self.intra.h_i = self.intra.h_i * self.intra.w / intra_w2

        # 4.3: Update intracellular volume, radius and surface area
        self.intra.w = intra_w2
        self.intra.radius = np.sqrt(self.intra.w / (self.intra.length * np.pi))
        self.intra.sa = 2 * np.pi * (self.intra.radius * self.intra.length)

        # 4.4: Update extracellular parameters
        if not self.infinite_bath:
            extra_w2 = self.extra.w - dw
            self.extra.na_i = self.extra.na_i * self.extra.w / extra_w2
            self.extra.k_i = self.extra.k_i * self.extra.w / extra_w2
            self.extra.cl_i = self.extra.cl_i * self.extra.w / extra_w2
            self.extra.x_i = self.extra.x_i * self.extra.w / extra_w2
            self.extra.w = extra_w2
            self.extra.radius = np.sqrt(self.extra.w / (self.extra.length * np.pi))
            self.extra.sa = 2 * np.pi * (self.extra.radius * self.extra.length)

        # Part 5: Raise exceptions

        if self.intra.cl_i < 0:
            print("Cl_i = " + str(self.intra.cl_i))
            print("d_Cl_i = " + str(self.intra.d_cl_i))
            raise Exception("[Cl-] < 0")

        if self.intra.k_i < 0:
            print("k_i = " + str(self.intra.k_i))
            print("d_k_i = " + str(self.intra.d_k_i))
            raise Exception("[K+] < 0 ")

    def first_save(self):

        # Open the new file and create the compartment and timing folders
        try:
            with h5py.File(self.new_file_name, mode='w') as self.hdf:
                self.hdf.create_group('COMPARTMENTS')
                group = self.hdf.get('COMPARTMENTS')

                group.create_group(name=self.intra.name)
                subgroup = group.get(self.intra.name)
                data_array = self.intra.get_array(self.run_t)
                subgroup.create_dataset(name=str(self.step_num), data=data_array)

                group.create_group(name=self.extra.name)
                subgroup = group.get(self.extra.name)
                data_array = self.extra.get_array(self.run_t)
                subgroup.create_dataset(name=str(self.step_num), data=data_array)

                self.hdf.create_group("TIMING")
                timing = self.hdf.get("TIMING")
                timing.create_dataset("DT", data=self.dt)
                timing.create_dataset("TOTAL_T", data=self.total_t)
                timing.create_dataset("INTERVALS", data=self.total_intervals)

                print("simulation file ('" + self.new_file_name + "') created in base directory")
        except OSError as e:
            raise Exception(f"Error creating file: {e}")

    def save_to_HDF(self):
        with h5py.File(self.new_file_name, mode='a') as self.hdf:
            group = self.hdf.get('COMPARTMENTS')

            subgroup = group.get(self.intra.name)
            data_array = self.intra.get_array(self.run_t)
            subgroup.create_dataset(name=str(self.step_num), data=data_array)

            subgroup = group.get(self.extra.name)
            data_array = self.extra.get_array(self.run_t)
            subgroup.create_dataset(name=str(self.step_num), data=data_array)

    def print_status(self, index):
        # index is the position in the status array

        sim_current_duration = time.time() - self.sim_start_t
        print(str(self.status_percentages[index] * 100) + " % complete in " + str(
            round(sim_current_duration, 2)) + " s")
        #some readout
        print("Vm: " + str(self.intra.v*1000) + "mV")
        print("Ecl: " + str(self.intra.E_cl*1000) + "mV")
        print("Ehco3-: " + str(self.intra.E_hco3*1000) + "mV")
        print("Eh: " + str(self.intra.E_h*1000) + "mV")
        print("Eh2 : " + str( -1 * RTF * np.log(self.intra.h_i / self.extra.h_i)*1000) + "mV")
        
        
        print("pH: " + str(-np.log10(self.intra.h_i)))
        print("Volume: " + str(self.intra.w))

        if index == 2:  # time to complete 1%
            hundred_percent_t = sim_current_duration * 100
            print("Estimated time to complete :" + str(
                round(hundred_percent_t / 60, 2)) + " minutes")

    def run_simulation(self):
        self.sim_start_t = time.time()
        self.first_save()

        while self.step_num < self.total_steps:

            self.sim_step()

            if self.step_num in self.status_steps:
                self.print_status(self.status_steps.index(self.step_num))

            if self.step_num == self.interval_arr[self.interval_num]:
                self.interval_num += 1
                if self.step_num != 0:
                    self.save_to_HDF()

            self.step_num = self.step_num + 1
            self.run_t = self.run_t + self.dt

        print("100.0 % complete in " + str(round(time.time() - self.sim_start_t, 2)) + " s")
        print("Simulation over")
        print("HDF file location: " + self.new_file_name)
