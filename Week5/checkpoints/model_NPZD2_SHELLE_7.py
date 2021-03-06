
'''
Fennel et al (2006) Nitrogen cycling in the Middle Atlantic Bight: Results from a 
three-dimensional model and implications for the North Atlantic nitrogen budget.
GLOBAL BIOGEOCHEMICAL CYCLES, VOL. 20, GB3007, doi:10.1029/2005GB002456

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
    # NPZD2
    par['mu0']   = 0.69  
    par['kNO3']  = 0.5    
    par['kNH4']  = 0.5  
    par['alpha'] = 0.125  
    par['gmax']  = 0.6
    par['kP']    = 0.44
    par['mP']    = 0.15    
    par['tau']   = 0.005 
    par['thetaMax'] = 0.053
    par['beta']  = 0.75 
    par['lBM']   = 0.1    
    par['lE']    = 0.1
    par['mZ']    = 0.25
    par['rSD']   = 0.3
    par['rLD']   = 0.1 
    par['nmax']  = 0.05
    par['kI']    = 0.1
    par['I0']    = 0.0095
    # SHELLE
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
    #NDPZD2
    InitCond = {}
    InitCond['Phy']  = 0.2
    InitCond['Zoo']  = 0.1
    InitCond['SDet'] = 1. 
    InitCond['LDet'] = 1.
    InitCond['NH4']  = 0.1
    InitCond['NO3']  = 7.
    InitCond['Temp'] = 6.
    # SHELLE
    InitCond['Soma'] = 0.01
    InitCond['Gonad'] = 0.
    InitCond['Salt'] = 30. #Salinity
    InitCond['Oxy'] = 30. #Oxygen
    return days, dt, par, InitCond
    



    
def run(days,dt,InitCond,par):
    '''
    This is your model. Do a brief description.

    '''
    # Import libraries
    import numpy as np
    
    # Setup the framework (calculate timestemps, create zero vectors, create time vector)
    NoSTEPS = int(days / dt) # Calculates the number of steps by dividing days by dt and rounding down
    time = np.linspace(0,days,NoSTEPS) # Makes and vector array of equally spaced numbers from zero to "days"
    
    # Create zero-vectors
    # NPZD2
    Phy = np.zeros((NoSTEPS,),float) # makes a vector array of zeros (size: NoSTEPS rows by ONE column)
    Zoo = np.zeros((NoSTEPS,),float) # same as above
    SDet = np.zeros((NoSTEPS,),float) # Biomass - same as above 
    LDet = np.zeros((NoSTEPS,),float) # same as above
    NH4 = np.zeros((NoSTEPS,),float) # same as above
    NO3 = np.zeros((NoSTEPS,),float) # same as above
    I  = np.zeros((NoSTEPS,),float) # same as above
    
    mu = np.zeros((NoSTEPS,),float) # same as above
    f_I = np.zeros((NoSTEPS,),float) # same as above
    L_NO3 = np.zeros((NoSTEPS,),float) # same as above
    L_NH4 = np.zeros((NoSTEPS,),float) # same as above
    TotN = np.zeros((NoSTEPS,),float) # same as above

    # SHELLE
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

    
    # Creating sunlight
    for i in range(len(I)):
        I[i] = 10 * np.sin((2*np.pi*time[i])/1) + \
               8 * np.sin((2*np.pi*time[i])/365)
        # We can't have negative light... so negatives are made zero 
        if I[i] < 0:
            I[i] = 0.0000001
    
    
    # Initializing with initial conditions
    # NPZD2
    Phy[0] = InitCond['Phy']
    Zoo[0] = InitCond['Zoo']
    SDet[0] = InitCond['SDet']
    LDet[0] = InitCond['LDet']
    NH4[0] = InitCond['NH4']
    NO3[0] = InitCond['NO3']
    
    # Initializing TotN
    TotN[0] = Phy[0] + Zoo[0] + SDet[0] + LDet[0] + NH4[0] + NO3[0]
    
    Temp = np.ones((NoSTEPS,),float) * InitCond['Temp'] #Temperature
    
    # SHELLE
    Soma[0] = InitCond['Soma']
    Gonad[0] = InitCond['Soma']
    B[0] = InitCond['Soma'] + InitCond['Gonad']
    Spawning[0] = 0.
    Salt = np.ones((NoSTEPS,),float) * InitCond['Salt']#Salinity
    Oxy = np.ones((NoSTEPS,),float) * InitCond['Oxy'] #Oxygen
    
    
    
    # *****************************************************************************
    # MAIN MODEL LOOP *************************************************************
    for t in range(0,NoSTEPS-1):
        # NPZD2
        muMax = par['mu0'] * (1.066 ** Temp[t]) # text
        
        f_I[t] = (par['alpha']*I[t])/(np.sqrt(muMax**2+((par['alpha']**2)*(I[t]**2)))) #Eq5
        
        L_NO3[t] = (NO3[t]/(par['kNO3']+NO3[t])) * (1/(1+(NH4[t]/par['kNH4']))) #Eq3
        
        L_NH4[t] = NH4[t]/(par['kNH4']+NH4[t]) # Eq 4
    
        mu[t] =muMax * f_I[t] * (L_NO3[t] + L_NH4[t]) # Eq2
    
        g = par['gmax'] * (Phy[t]**2/(par['kP']+(Phy[t]**2)))
        
        n = par['nmax'] * (1 - max(0,(I[t]-par['I0'])/(par['kI']+I[t]-par['I0'])))
    
        dPhydt = (mu[t] * Phy[t]) - \
                 (g  * Zoo[t]) - \
                 (par['mP'] * Phy[t]) - \
                 (par['tau']*(SDet[t]+Phy[t])*Phy[t]) # Eq1
                 
        dZoodt = (g * par['beta'] * Zoo[t]) - \
                 (par['lBM']*Zoo[t]) - \
                 (par['lE']*((Phy[t]**2)/(par['kP']+(Phy[t]**2)))*par['beta']*Zoo[t]) - \
                 (par['mZ']*(Zoo[t]**2))#Eq10
                 
        dSDetdt = (g * (1-par['beta']) * Zoo[t]) + \
                  (par['mZ']*(Zoo[t]**2)) + \
                  (par['mP'] * Phy[t]) - \
                  (par['tau']*(SDet[t]+Phy[t])*SDet[t]) - \
                  (par['rSD']*SDet[t])
                  
        dLDetdt = (par['tau']*((SDet[t]+Phy[t])**2)) - \
                  (par['rLD']*LDet[t])
                  
        dNO3dt = -(muMax * f_I[t] * L_NO3[t] * Phy[t]) + \
                  (n * NH4[t])
                 
        dNH4dt = -(muMax * f_I[t] * L_NH4[t] * Phy[t]) - \
                  (n * NH4[t]) + \
                  (par['lBM'] * Zoo[t]) + \
                  (par['lE']*((Phy[t]**2)/(par['kP']+(Phy[t]**2)))*par['beta']*Zoo[t]) + \
                  (par['rSD']*SDet[t]) + \
                  (par['rLD']*LDet[t])
                  
        # SHELLE
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
        # NPZD2
        Phy[t+1]  = Phy[t]  + (dPhydt * dt)
        Zoo[t+1]  = Zoo[t]  + (dZoodt * dt)
        SDet[t+1] = SDet[t] + (dSDetdt * dt)
        LDet[t+1] = LDet[t] + (dLDetdt * dt)
        NH4[t+1]  = NH4[t]  + (dNH4dt * dt)
        NO3[t+1]  = NO3[t]  + (dNO3dt * dt)
        if NO3[t+1] <= 0.0001:
            offset = NO3[t+1]
            NH4[t+1] = NH4[t+1] + offset
            NO3[t+1] = NO3[t+1] - offset
        
        # SHELLE
        Soma[t+1] = Soma[t] + (dSomadt * dt)
        Gonad[t+1] = Gonad[t] + (dGonaddt * dt)  - Spawning[t]
        B[t+1] = Soma[t+1] + Gonad[t+1]
            
        # Estimate Total Nitrogen
        TotN[t+1] = Phy[t+1] + Zoo[t+1] + SDet[t+1] + LDet[t+1] + NH4[t+1] + NO3[t+1]
    # end of main model LOOP*******************************************************
    # *****************************************************************************

    # Pack output into dictionary
    output = {}
    output['time'] = time
    # NPZD2
    output['Phy'] = Phy
    output['Zoo'] = Zoo
    output['SDet'] = SDet
    output['LDet'] = LDet
    output['NH4'] = NH4
    output['NO3'] = NO3
    output['mu'] = mu
    output['f_I'] = f_I
    output['L_NO3'] = L_NO3
    output['L_NH4'] = L_NH4
    output['TotN'] = TotN
    #SHELLE
    output['Soma'] = Soma
    output['Gonad'] = Gonad
    output['B'] = B
    output['Spawning'] = Spawning
    output['Temp'] = Temp
    output['Salt'] = Salt

    print('Model run: DONE!!!')
    return output


    
def plot(output):
    '''
    Script to make plots
    '''
    # Import libraries
    import matplotlib.pyplot as plt
    
    # Plotting
    fig, (ax, ax2) = plt.subplots(2,1,figsize=(15, 8))

    ax.plot(output['time']/365,output['f_I'],'r-')
    ax.plot(output['time']/365,output['mu'],'g-')
    ax.plot(output['time']/365,output['L_NO3'],'b-')
    ax.plot(output['time']/365,output['L_NH4'],'k-')
    ax.set_xlabel('Time (years)')
    ax.set_ylabel('Limitations (dimensionless)')
    ax.set_title('Fennel et al 2006 Model')
    ax.legend(['f_I','Mu','L_NO3','L_NH4'])
    
    ax2.plot(output['time']/365,output['Phy'],'g-')
    ax2.plot(output['time']/365,output['Zoo'],'r-')
    ax2.plot(output['time']/365,output['SDet'],'k-')
    ax2.plot(output['time']/365,output['LDet'],'k-.')
    ax2.plot(output['time']/365,output['NH4'],'m-')
    ax2.plot(output['time']/365,output['NO3'],'c-')
    ax2.plot(output['time']/365,output['TotN'],'y-')
    ax2.set_xlabel('Time (years)')
    ax2.set_ylabel('Nitrogen (mmol N m$^{-3}$)')
    plt.legend(['Phy','Zoo','SDet','LDet','NH4','NO3','TotN'])
    plt.show()
    return
    
if __name__ == '__main__':
    days, dt, par, InitCond = load_defaults()
    output = run(days,dt,InitCond,par)
    plot(output)