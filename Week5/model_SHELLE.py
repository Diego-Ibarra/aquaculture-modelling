
'''
Here you can put a brief description of the model
and the reference of the paper you got it

Ibarra, Fennel, Cullen (2014) Coupling 3-D Eulerian bio-physics (ROMS) with individual-based
shellfish ecophysiology (SHELL-E): A hybrid model for carrying
capacity and environmental impacts of bivalve aquaculture. Ecological Modelling 273 (2014) 63- 78

'''


def load_defaults():
    '''
    This function creates a dictionaries called "par" and "InitCond"
    and pre-loads them with all the default 
    parameters and initial conditions, respectively.
    Also outputs days and dt
    '''
    # Framework
    days = 365 * 3 # Three year
    dt   = 0.01 # units: days    
    
    # Parameters
    par = {}
    par['AE_P']    = 0.9  
    par['AE_D']    = 0.2    
    par['AE_Z']    = 0.3  
    par['Bpub']    = 0.43  
    par['Fmax_ref'] = 0.025
    par['GT']       = 0.44
    par['KTempH']   = 0.1    
    par['KTempL']   = 0.5 
    par['KSaltL']   = 0.25
    par['KOxyL']    = 0.02 
    par['KFood']    = 1.    
    par['KRE']   = 0.86
    par['OxyL']  = 17.5
    par['Rm']    = 0.002
    par['SaltL'] = 10.
    par['TempH'] = 25.
    par['TempL'] = -4.
    par['beta']  = 0.12
    par['epsilonP'] = 1.
    par['epsilonD'] = 0.5
    par['epsilonZ'] = 0.3
    
    # Initial conditions
    InitCond = {}
    InitCond['Soma'] = 0.01
    InitCond['Gonad'] = 0.
    InitCond['Temp'] = 10. #Temperature
    InitCond['Salt'] = 30. #Salinity
    InitCond['Oxy'] = 30. #Oxygen
    InitCond['Phy'] = 2 #Phyto
    InitCond['Zoo'] = 0.5 #Zoo
    InitCond['SDet'] = 0.5 #SDet
    
    return days, dt, par, InitCond
    



    
def run(days,dt,InitCond,par):
    '''
    This is your model. Do a brief description.
    
    INPUTS:
        days: number of days of simulation
        dt: time steps (units: days)
        InitCond: Dictionary with all initial conditions
        par: Dictionary with all model parameters
        
    OUTPUTS:
        var1: name (units)
        var2: name (units)
        var3: name (units)
    
    Don't forget to reference the paper where you got it    
    '''
    # Import libraries
    import numpy as np
    
    # Setup the framework (calculate timestemps, create zero vectors, create time vector)
    NoSTEPS = int(days / dt) # Calculates the number of steps by dividing days by dt and rounding down
    time = np.linspace(0,days,NoSTEPS) # Makes and vector array of equally spaced numbers from zero to "days"
    
    # Create zero-vectors
    Soma = np.zeros((NoSTEPS,),float) # makes a vector array of zeros (size: NoSTEPS rows by ONE column)
    Gonad = np.zeros((NoSTEPS,),float) # same as above
    B = np.zeros((NoSTEPS,),float) # Biomass - same as above 
    L_Temp = np.zeros((NoSTEPS,),float) # same as above
    L_Salt = np.zeros((NoSTEPS,),float) # same as above
    L_Oxy = np.zeros((NoSTEPS,),float) # same as above
    L_Food = np.zeros((NoSTEPS,),float) # same as above
    F = np.zeros((NoSTEPS,),float) # same as above
    A = np.zeros((NoSTEPS,),float) # same as above
    R = np.zeros((NoSTEPS,),float) # same as above
    RE = np.zeros((NoSTEPS,),float) # same as above
    Spawning = np.zeros((NoSTEPS,),float) # same as above
    
    # Initializing with initial conditions
    Soma[0] = InitCond['Soma']
    Gonad[0] = InitCond['Soma']
    B[0] = InitCond['Soma'] + InitCond['Gonad']
    Spawning[0] = 0.
    Temp = np.ones((NoSTEPS,),float) * InitCond['Temp'] #Temperature
    Salt = np.ones((NoSTEPS,),float) * InitCond['Salt']#Salinity
    Oxy = np.ones((NoSTEPS,),float) * InitCond['Oxy'] #Oxygen
    Phy = np.ones((NoSTEPS,),float) * InitCond['Phy'] #Phyto
    Zoo = np.ones((NoSTEPS,),float) * InitCond['Zoo'] #Zoo
    SDet = np.ones((NoSTEPS,),float) * InitCond['SDet'] #SDet
    
    
    
    # *****************************************************************************
    # MAIN MODEL LOOP *************************************************************
    for t in range(0,NoSTEPS-1):
    # Calculate Temperature Limitation
        L_Temp[t] = min(max(0.,1.-np.exp(-par['KTempL']*(Temp[t]-par['TempL']))), \
                     max(0.,1.+((1.-np.exp(par['KTempH']*Temp[t]))/(np.exp(par['KTempH']*par['TempH'])-1.))))
        
        # Calculate Salinity Limitation
        L_Salt[t] = max(0.,1.-np.exp(-par['KSaltL']*(Salt[t]-par['SaltL'])))
        
        # Calculate Oxygen Limitation
        L_Oxy[t] = max(0.,1.-np.exp(-par['KOxyL']*(Oxy[t]-par['OxyL'])))
        
        # Calculate Oxygen Limitation
        L_Food[t] = (Phy[t]+Zoo[t]+SDet[t])/(par['KFood']+Phy[t]+Zoo[t]+SDet[t])
        
        # Calculate Filtration rate
        Fmax  = par['Fmax_ref']*(B[t]**(2./3.))
        
        F[t] = Fmax * L_Temp[t] * L_Salt[t] * L_Oxy[t] * L_Food[t]
        
        A[t] = F[t] * ((par['epsilonP']*par['AE_P']*Phy[t])+ \
                       (par['epsilonZ']*par['AE_Z']*Zoo[t])+ \
                       (par['epsilonD']*par['AE_D']*SDet[t]))
        
        R[t] = (par['Rm']*B[t]) + (par['beta']*A[t])
        
        RE[t] = max(0., (B[t]-par['Bpub'])/(par['KRE'] + B[t] - (2.*par['Bpub'])))
        
        # Spawning
        if Gonad[t]/B[t] < par['GT']:
            Spawning[t] = 0.
            dGonaddt = (A[t]-R[t]) * RE[t]
            dSomadt =  (A[t]-R[t]) * (1.-RE[t])
        elif Gonad[t]/B[t] >= par['GT']:         
            Spawning[t] = Gonad[t]
            dGonaddt = 0.
            dSomadt = A[t]-R[t]
        else:
            dGonaddt = 0.
            dSomadt = 0.   
        
        # Update and step ----------------------------------------------------
        Soma[t+1] = Soma[t] + (dSomadt * dt)
        Gonad[t+1] = Gonad[t] + (dGonaddt * dt)  - Spawning[t]
        B[t+1] = Soma[t+1] + Gonad[t+1]
    # end of main model LOOP*******************************************************
    # *****************************************************************************

    # Pack output into dictionary
    output = {}
    output['time'] = time
    output['Soma'] = Soma
    output['Gonad'] = Gonad
    output['B'] = B
    output['Spawning'] = Spawning
    output['Temp'] = Temp
    output['Salt'] = Salt
    output['Oxy'] = Oxy
    output['Phy'] = Phy
    output['Zoo'] = Zoo
    output['SDet'] = SDet

    print('Model run: DONE!!!')
    return output


    
def plot(output):
    '''
    Script to make plots
    '''
    # Import libraries
    import matplotlib.pyplot as plt
    
    # Plotting
    fig, (ax) = plt.subplots(1,1,figsize=(15, 4))
    ax.plot(output['time']/365,output['Soma'],'b-')
    ax.plot(output['time']/365,output['Gonad'],'g.')
    ax.plot(output['time']/365,output['B'],'r-')
    ax.set_xlabel('Time (years)')
    ax.set_ylabel('Nitrogen (mmol N m$^{-3}$)')
    ax.set_title('Mussel_IbarraEtal2014 Model Simulation')
    plt.legend(['Soma','Gonad','Biomass'],loc='upper left')
    plt.show()
    return
    
if __name__ == '__main__':
    days, dt, par, InitCond = load_defaults()
    output = run(days,dt,InitCond,par)
    plot(output)