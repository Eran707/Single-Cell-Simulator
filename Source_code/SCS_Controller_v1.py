# ---------------------------------------------------
####################################################
# 20 November 2023
# Contributors: Eran Shorer, Kira Dusterwald, Christopher Currin, Joseph Raimondo

# SCS - SINGLE CELL SIMULATOR

# This software offers the ability to simulate neural dynamics
# The script is the controller which allows the user to setup and run the simulation
# The SCS_Simulator script provides the backend functionality

# The SCS saves the results of the simulation in an HDF file
# When starting the simulation one can use the settings of a previous simulation or start a simulation from scratch

#####################################################
# ---------------------------------------------------

from SCS_Simulator_v1 import Simulator

# STEP 1: ESTABLISHING SIMULATION TYPE

sim_type = "New"  # Simulation type either set as "New" or "Extend"

## Baseline

new_file_name = "default_oldVersion"
old_file_name = "rate_constants_SS_default_z_ratio_constants_exch"  # Only needed when the sim_type is Extend

## zChange:

# new_file_name = "z_change_oldVersion"

## KCC2:

#new_file_name = "KCC2_change_oldVersion"

# STEP 2: SET SIMULATION TIMING


total_t = 6000  # total simulation time in seconds
dt = 1e-5  # simulation time step in seconds
intervals = 600000  # number of times the results of the simulation will be saved to the HDF file

# STEP 3: SET SIMULATION SETTINGS

infinite_bath = True
dynamic_ATPase = True

# STEP 4: SET ION CONCENTRATIONS (ONLY FOR NEW SIMULATIONS)


na_i = 15e-3
k_i = 122.6e-3
cl_i = 5.3e-3
hco3_i = 9.7e-3
# h_i = 6.31e-8 # we did not have this in the older version
X_i = 144.4e-3
Total_intra_charge = na_i + k_i - (cl_i + hco3_i + X_i * 0.85)
intracellular_params = {"name": "1_Intracellular", "radius": 5e-5, "length": 25e-5,
                        "na": na_i, "k": k_i, "cl": cl_i, "hco3": hco3_i, "X": X_i, "z": -0.85}
na_o = 145e-3
k_o = 3.5e-3
cl_o = 110e-3
hco3_o = 31e-3
# h_o = 3.98e-8 - h_o
X_o = 7.5e-3
Total_extra_charge = na_o + k_o - (cl_o + hco3_o + X_o * 1)
extracellular_params = {"name": "2_Extracellular", "radius": 5e-5, "length": 25e-5,
                        "na": na_o, "k": k_o, "cl": cl_o, "hco3": hco3_o, "X": 7.5e-3, "z": -1}

# STEP 5: ESTABLISH BASIC SIMULATION PARAMETERS

if sim_type == "New":
    sim = Simulator(sim_type="New", new_file_name=new_file_name)
    sim.set_intracellular_properties(intracellular_params)
    sim.set_extracellular_properties(extracellular_params)

elif sim_type == "Extend":
    sim = Simulator(sim_type="Extend", new_file_name=new_file_name, old_file_name=old_file_name)

else:
    print("sim_type needs to either be 'New' or 'Extend'")

sim.set_timing(total_t=total_t, dt=dt, total_intervals=intervals)
sim.infinite_bath = infinite_bath
sim.dynamic_ATPase = dynamic_ATPase
# sim.set_kcc2_off()


# STEP 6: ADD ADDITIONAL COMPONENTS TO THE SIMULATION

## Options available: Z_change; X_change; Add_synapse; Add_current; Change_KCC2


## Baseline

#sim.add_synapse(start_t=2000, tau=50e-3, max_g=1e-8)

## zChange

# sim.add_synapse(start_t=2000, tau=50e-3, max_g=1e-8)
# sim.set_z_change(start_t=2500,end_t=3000,z_change_amount=-0.2,adjust_cl=False)
# sim.add_synapse(start_t=3500, tau=50e-3, max_g=1e-8)

## KCC2

#sim.add_synapse(start_t=2000, tau=50e-3, max_g=1e-8)
#sim.set_KCC2_change(start_t=2500, end_t=3000, final_KCC2_value=0)
#ssim.add_synapse(start_t=5000, tau=50e-3, max_g=1e-8)

sim.run_simulation()
