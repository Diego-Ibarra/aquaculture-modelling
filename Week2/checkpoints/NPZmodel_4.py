# Import libraries
import numpy as np

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

#%% Chores (calculate timestemps, create zero vectors, create time vector)
NoSTEPS = int(days / dt) # Calculates the number of steps by dividing days by dt and rounding down
time = np.linspace(0,days,NoSTEPS) # Makes and vector array of equally spaced numbers from zero to "days"
N = np.zeros((NoSTEPS,),float) # makes a vector array of zeros (size: NoSTEPS rows by ONE column)
P = np.zeros((NoSTEPS,),float) # same as above
Z = np.zeros((NoSTEPS,),float) # same as above

#%% Initializing with initial conditions
P[0] = Pinitial
N[0] = Ninitial
Z[0] = Zinitial

#%% *****************************************************************************
# MAIN MODEL LOOP *************************************************************
for t in range(0,NoSTEPS-1):
    L_N = N[t]/(K+N[t]) # Calculate Limitation due to (low) nutrients on pythoplankton
    L_P = 1-np.exp(-Lambda*P[t]) # Calculate Limitation due to (low) pythoplankton on zooplankton
    
    # Estimate model state at time t+1 
    P[t+1] = P[t] + (((mu*L_N*P[t]) - (gamma*L_P*Z[t])- (mP*P[t])) * dt)
    N[t+1] = N[t] + (((mP*P[t]) + (mZ*Z[t]) - (mu*L_N*P[t])) * dt)
    Z[t+1] = Z[t] + (((gamma*L_P*Z[t]) - (mZ*Z[t])) * dt)
# end of main model LOOP*******************************************************
# *****************************************************************************