#%% Framework
days = 10
dt   = 0.01 # units: days

# Parameters
mu     = 0.5  # growth rate of pythoplankton (units: d^-1)
K      = 1    # half-saturation for nutrient absoption by pythoplankton (units: mmolN m^-3)
gamma  = 0.4  # growth rate of zooplankton (units: d^-1)
Lambda = 0.1  # initial slope of ingestion saturation of zooplankton (units: mmol N^-1 m^3)
mP     = 0.07 # mortality rate of pythoplankton (units: d^-1)
mZ     = 0.03 # mortality rate of zooplankton (units: d^-1)

# Initial conditions
Pinitial = 4  # initial pythoplankton concentration (units: mmolN m^-3)
Ninitial = 10 # initial nutrients concentration (units: mmolN m^-3)
Zinitial = 2  # initial zooplankton concentration (units: mmolN m^-3)
