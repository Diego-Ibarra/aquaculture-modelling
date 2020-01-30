def load_defaults():
    print('Loading defaults...')
    # Framework
    days = 365 * 3 # One year
    dt   = 0.01 # units: days    
    
    # Parameters
    par = {}
    par['AE_P']  = 0.9  
    par['AE_D']  = 0.2 
    par['AE_Z']  = 0.3
    par['Fmax_ref']= 0.025 
    par['epsilonP'] = 1. 
    par['epsilonD'] = 0.5 
    par['epsilonZ'] = 0.3
    par['KTempH']= 0.1  
    par['KTempL']= 0.5
    par['TempH'] = 25. 
    par['TempL'] = -4. 
    par['KSaltL']= 0.25 
    par['SaltL'] = 10.
    par['OxyL']  = 17.5 
    par['KOxyL'] = 0.02
    par['KFood'] = 1.  
    par['Rm']    = 0.002 
    par['beta']  = 0.12
    par['Bpub']    = 0.43
    par['KRE']   = 0.86
    par['GT']    = 0.44
    
    # Initial conditions
    InitCond = {}
    InitCond['Soma'] = 0.01 
    InitCond['Gonad'] = 0.0
    InitCond['Phy'] = 0.8 
    InitCond['Zoo'] = 0.3 
    InitCond['SDet'] = 0.2
    InitCond['Temp'] = 10 
    InitCond['Salt'] = 30 
    InitCond['Oxy'] = 340.
    return  days, dt, par, InitCond
    
def run(days, dt, par, InitCond):
    print('Running model...')
    # Import libraries
    import numpy as np
    
    # Setup the framework 
    NoSTEPS = int(days / dt) # Calculates the number of steps 
    time = np.linspace(0,days,NoSTEPS) # Makes vector array of equally spaced numbers 
    
    # Create arrays of zeros
    B = np.zeros((NoSTEPS,),float) # Biomass 
    Soma = np.zeros((NoSTEPS,),float) 
    Gonad = np.zeros((NoSTEPS,),float)
    
    # Initializing with initial conditions
    Soma[0] = InitCond['Soma'] 
    Gonad[0] = InitCond['Gonad'] 
    B[0] = InitCond['Soma'] + InitCond['Gonad'] 
    
    # *****************************************************************************
    # MAIN MODEL LOOP *************************************************************
    for t in range(0,NoSTEPS-1):
        
        Fmax = par['Fmax_ref']*(B[t]**(2./3.)) # Eq.5
        
        # Eq.6 - Temperature Limitation  
        L_Temp = min(max(0.,1.-np.exp(-par['KTempL']*(InitCond['Temp']-par['TempL']))), \
                     max(0.,1.+((1.-np.exp(par['KTempH']*InitCond['Temp']))/(np.exp(par['KTempH']*par['TempH'])-1.))))
        
        # Eq.7 - Salinity Limitation 
        L_Salt = max(0.,1.-np.exp(-par['KSaltL']*(InitCond['Salt']-par['SaltL'])))
        
        # Eq.8 - Oxygen Limitation
        L_Oxy = max(0.,1.-np.exp(-par['KOxyL']*(InitCond['Oxy']-par['OxyL'])))
        
         # Eq.9 - Food Limitation
        L_Food = (InitCond['Phy']+InitCond['Zoo']+InitCond['SDet'])/ \
                 (par['KFood']+InitCond['Phy']+InitCond['Zoo']+InitCond['SDet']) 
        
        F = Fmax * L_Temp * L_Salt * L_Oxy * L_Food #Eq 4 

        # Eq.3 # 
        A = F * ((par['epsilonP']*par['AE_P']*InitCond['Phy'])+ \
                    (par['epsilonZ']*par['AE_Z']*InitCond['Zoo'])+ \
                    (par['epsilonD']*par['AE_D']*InitCond['SDet']))
        
        R = (par['Rm']*B[t]) + (par['beta']*A)
        
        # Dynamic mass allocation - Eq. 13 
        RE = max(0., (B[t]-par['Bpub'])/(par['KRE'] + B[t] - (2.*par['Bpub']))) 
        
        
        # Spawning Eq. 14 
        if Gonad[t]/B[t] < par['GT']:
            Spawning = 0.
        elif Gonad[t]/B[t] >= par['GT']:
            Spawning = Gonad[t]
            Gonad[t] = 0.
            
        
        #dBdt = (A - R) - Spawning # Eq.2 
        dSomadt = (A-R) * (1.-RE) # Eq. 11
        
        # Note that I added a "max" function to prevent Gonad to go negative 
        dGonaddt = max(0.,((A-R) * RE) - Spawning) #Eq 12. 
        
        # Update and step ------------------------------
        Soma[t+1] = Soma[t] + (dSomadt * dt) 
        Gonad[t+1] = Gonad[t] + (dGonaddt * dt) 
        B[t+1] = Soma[t+1] + Gonad[t+1] 
        
    # end of main model LOOP*******************************************************
    # *****************************************************************************

    # Pack output into dictionary
    output = {}
    output['time'] = time
    output['B'] = B 
    output['Soma'] = Soma 
    output['Gonad'] = Gonad 
    output['L_Temp'] = L_Temp 
    output['L_Salt'] = L_Salt 
    output['L_Oxy'] = L_Oxy
    output['L_Food'] = L_Food 
    
    print('Model run: DONE!!!')
    return  output

def plot(output):
    import matplotlib.pyplot as plt 
    # Plotting                      
    fig, (ax) = plt.subplots(1,1)
    ax.plot(output['time']/365,output['B'],'r-') 
    ax.plot(output['time']/365,output['Soma'],'b-') 
    ax.plot(output['time']/365,output['Gonad'],'g.')
    ax.legend(['B', 'Soma', 'Gonad']) 
    ax.set_ylabel('Nitrogen (mmol N)') 
    ax.set_xlabel('Time (years)')
    plt.show()          
    print('L_Temp = ' + str(output['L_Temp']))
    print('L_Salt = ' + str(output['L_Salt']))
    print('L_Oxy  = ' + str(output['L_Oxy']))
    print('L_Food  = ' + str(output['L_Food']))
    return

if __name__ == "__main__":
    print('Executing my_module.py')
    print('--------------------')
    
    days, dt, par, InitCond = load_defaults()
    output = run(days, dt, par, InitCond)
    plot(output)
    
    print('--------------------')
    