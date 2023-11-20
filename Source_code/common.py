"""
Simulation parameters that are common to all simulations

Created and modified by: K Dusterwald; C Currin; EF Shorer
"""
###########################
# Constants
###########################


R = 8.31446
F: float = 96485.33  # Faraday's constant        C mol-1
k = 1.38e-23  # Boltzmann constant        J K-1
q = 1.602176620898e-19  # elementary charge         C
Na = 6.022e23  # Avogadro's constant       mol-1

T = 37 + 273.15
RTF = R * T / F
RT = R * T

vw = 0.018  # partial molar volume of water, dm3/mol
pw = 0.0015  # osmotic permeability, biological membrane (muscle? unknown), dm s
km = 5 * 10 ** (-14)  # extensional rigidity of RBC at 23 deg, Mohandas and Evans (1994), N/dm

##########################
# Cellular environment
##########################

# Membrane capacitance
cm = 2e-4  # (F/dm^2)


# Ion valency
val = {"na": 1, "k": 1, "cl": -1, "x": -0.85}

# Electrodiffusion coefficients
diff_constants = {"na": 1.33e-7, "k": 1.96e-7, "cl": 2.03e-7, "x": 0}  # dm2/s

# Acid base

pKa = 6.1 # % Acid dissociation constant for H2CO3
kH = 3.1e-2 # % Henry's law constant for CO2 in water (M/atm)
P_CO2 = 1 * 5/100 # Calculate PCO2 from given percent CO2Convert % CO2 to atmospheres
h2co3_i = kH * P_CO2 # Calculate H2CO3 concentration using Henry's law
h_i = 10**(-7.2) # [H+] in mM based on a pH of 7.2

kf = 10 **3 #  forward rate constant for carbonic acid equation


###########################
# Ion channels and pumps
###########################

# Leak conductances:
#gna = (2e-3) / (F) #default
#gna = (3.1e-3) / (F)
gna = (2.5e-3) / (F)
gk = (7e-3) / (F)
gcl = (2e-3) / (F)  # gna,gk,gcl: conductances in mS/cm^2 conv to S/dm^2 (10^-3/10^-2) - corrected for neuron
ghco3 = gcl*0.2 # ghco3 is 20% of gcl

# KCC2 and ATPase pumps
p_kcc2 = (2e-3) / (F) #default
p_atpase = 0.1 / F  # C/(dm2·s)


