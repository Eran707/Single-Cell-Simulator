# ---------------------------------------------------
####################################################
# 30 October 2023
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

sim_type = "Extend"    # Simulation type either set as "New" or "Extend"

new_file_name = "SS_default_synapse"
old_file_name = "SS_default"  # Only needed when the sim_type is Extend

# STEP 2: SET SIMULATION TIMING

total_t = 10      # total simulation time in seconds
dt = 1e-5          # simulation time step in seconds
intervals = 10000   # number of times the results of the simulation will be saved to the HDF file

# STEP 3: SET SIMULATION SETTINGS

infinite_bath = True
dynamic_ATPase = True


# STEP 4: SET ION CONCENTRATIONS (ONLY FOR NEW SIMULATIONS)

intracellular_params = {"name": "1_Intracellular", "radius": 5e-5,"length": 25e-5,
                        "na": 14e-3, "k":122.9e-3 , "cl":5.2e-3, "hco3":15e-3, "X":154.9e-3,"z":-0.85 }

extracellular_params = {"name": "2_Extracellular", "radius": 5e-5, "length": 25e-5,
                        "na": 145e-3,"k":3.5e-3 , "cl":110e-3, "hco3":31e-3, "X":7.5e-3,  "z":-1}

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

#sim.set_kcc2_off()


# STEP 6: ADD ADDITIONAL COMPONENTS TO THE SIMULATION
# Options available: Z_change; X_change; Add_synapse; Add_current; Change_KCC2

sim.add_synapse(start_t=5, tau=250e-3, max_g=1e-9)

#sim.set_z_change(start_t=1,end_t=2,z_change_amount=-0.2)
#sim.set_KCC2_change(start_t=3000, end_t=6000, final_KCC2_value=0)

sim.run_simulation()



