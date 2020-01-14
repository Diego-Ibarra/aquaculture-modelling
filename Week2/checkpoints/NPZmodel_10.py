# Import libraries
import numpy as np
import matplotlib.pyplot as plt

#%% Framework
days = 365 * 2
dt   = 0.01 # units: days

# Parameters
mu     = 0.5  # growth rate of pythoplankton (units: d^-1)
K      = 1    # half-saturation for nutrient absoption by pythoplankton (units: mmolN m^-3)
gamma  = 0.4  # growth rate of zooplankton (units: d^-1)
Lambda = 0.1  # initial slope of ingestion saturation of zooplankton (units: mmol N^-1 m^3)
alpha  = 0.004# initial slope of P vs I curve (units: {mmolN m^-3}^-1)
psi    = 0.01    # half-saturation for phytoplankton absoption by mussels (units: mmolN m^-3)
mP     = 0.07 # mortality rate of pythoplankton (units: d^-1)
mZ     = 0.03 # mortality rate of zooplankton (units: d^-1)
mM     = 0.0001 # mortality rate of zooplankton (units: d^-1)
Km     = 0.1    # half-saturation for phytoplankton feeding by mussels due to mussel size (units: mmolN m^-3)
Kmp     = 3    # half-saturation for phytoplankton feeding by mussels (units: mmolN m^-3)

# Initial conditions
Pinitial = 4  # initial pythoplankton concentration (units: mmolN m^-3)
Ninitial = 10 # initial nutrients concentration (units: mmolN m^-3)
Zinitial = 2  # initial zooplankton concentration (units: mmolN m^-3)
Minitial = 10  # initial mussel concentration (units: mmolN m^-3)

#%% Chores (calculate timestemps, create zero vectors, create time vector)
NoSTEPS = int(days / dt) # Calculates the number of steps by dividing days by dt and rounding down
time = np.linspace(0,days,NoSTEPS) # Makes and vector array of equally spaced numbers from zero to "days"
N = np.zeros((NoSTEPS,),float) # makes a vector array of zeros (size: NoSTEPS rows by ONE column)
P = np.zeros((NoSTEPS,),float) # same as above
Z = np.zeros((NoSTEPS,),float) # same as above
M = np.zeros((NoSTEPS,),float) # same as above
TotNit = np.zeros((NoSTEPS,),float) # same as above
L_N = np.zeros((NoSTEPS,),float) # same as above
L_P = np.zeros((NoSTEPS,),float) # same as above
I  = np.zeros((NoSTEPS,),float) # same as above
L_I = np.zeros((NoSTEPS,),float) # same as above
L_M = np.zeros((NoSTEPS,),float) # same as above
L_MP = np.zeros((NoSTEPS,),float) # same as above

#%% Creating sunlight
for i in range(len(I)):
    I[i] = 600 * np.sin((2*np.pi*time[i])/1) + \
           500 * np.sin((2*np.pi*time[i])/365)
    
    # We can't have negative light... so negatives are made zero 
    if I[i] < 0:
        I[i] = 0

#%% Initializing with initial conditions
P[0] = Pinitial
N[0] = Ninitial
Z[0] = Zinitial
M[0] = Minitial


#%% *****************************************************************************
# MAIN MODEL LOOP *************************************************************
for t in range(0,NoSTEPS-1):
    L_N[t] = N[t]/(K+N[t]) # Calculate Limitation due to (low) nutrients on pythoplankton
    L_P[t] = 1-np.exp(-Lambda*P[t]) # Calculate Limitation due to (low) pythoplankton on zooplankton
    L_I[t] = 1-np.exp(-alpha*I[t])
    L_M[t] = M[t]/(Km+M[t]) # Calculate Limitation due to "old age" on mussels
    L_MP[t] = P[t]/(Kmp+P[t]) # Calculate Limitation due to "old age" on mussels
    
    # Estimate model state at time t+1 
    P[t+1] = P[t] + (((mu*L_N[t]*L_I[t]*P[t]) - (gamma*L_P[t]*Z[t])- (psi*L_M[t]*L_MP[t]*M[t]) - (mP*P[t])) * dt)
    N[t+1] = N[t] + (((mP*P[t]) + (mZ*Z[t]) + (mM*M[t]) - (mu*L_N[t]*L_I[t]*P[t])) * dt)
    Z[t+1] = Z[t] + (((gamma*L_P[t]*Z[t]) - (mZ*Z[t])) * dt)
    M[t+1] = M[t] + (((psi*L_M[t]*L_MP[t]*M[t]) - (mM*M[t])) * dt)
    TotNit[t+1] = P[t+1] + N[t+1] + Z[t+1] + M[t+1]
# end of main model LOOP*******************************************************
# *****************************************************************************
    
#%% Plotting
fig, (ax) = plt.subplots(1,1)
ax.plot(time,N,'b-')
ax.plot(time,P,'g-')
ax.plot(time,Z,'r-')
ax.plot(time,M,'k-')
ax.plot(time,TotNit,'c-')
ax.set_xlabel('Time (days)')
ax.set_ylabel('Nitrogen (mmol N m$^{-3}$)')
ax.set_title('NPZM Model Simulation')
plt.legend(['N','P','Z','M','TotNit'])
plt.show()

fig2, (ax) = plt.subplots(1,1)
ax.plot(time,L_I,'c-')
ax.plot(time,L_N,'r-')
ax.plot(time,L_P,'b-')
ax.plot(time,L_M,'k-')
ax.plot(time,L_MP,'m-')
ax.set_xlabel('Time (days)')
ax.set_ylabel('Limitation (dimensionless)')
ax.set_title('Time evolution of Limitations')
plt.legend(['L_N: Lim of N on Phy','L_P: Lim of P on Zoo','L_I: Lim of I on Phy','L_M: Lim of M on M','L_MP: Lim of P on M'])
plt.show()