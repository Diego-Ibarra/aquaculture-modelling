
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
    # Physical characteristics of embayment
    par['X'] = 1. # Basin length
    par['Y'] = 1. # Basin width
    par['Z'] = 1. # Basin depth
    par['V'] = par['X'] * par['Y'] * par['Z']
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
    par['betaSHELLE']  = 0.12
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
    InitCond['n_muss'] = 1.
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
    n_muss = np.zeros((NoSTEPS,),float) # same as above
    B_total = np.zeros((NoSTEPS,),float) # Total Biomass of mussels in Embayment 

    
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

    Temp = np.ones((NoSTEPS,),float) * InitCond['Temp'] #Temperature 
    
    # SHELLE
    Soma[0] = InitCond['Soma']
    Gonad[0] = InitCond['Soma']
    B[0] = InitCond['Soma'] + InitCond['Gonad']
    Spawning[0] = 0.
    Salt = np.ones((NoSTEPS,),float) * InitCond['Salt']#Salinity
    Oxy = np.ones((NoSTEPS,),float) * InitCond['Oxy'] #Oxygen
    n_muss[0] = InitCond['n_muss']
    B_total[0] = B[0] * n_muss[0] / par['V']

    
    
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
        
        R[t] = (par['Rm']*B[t]) + (par['betaSHELLE']*A[t])
        
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
            Zoo[t] =  Zoo[t] + ((Spawning[t]*n_muss[t])/par['V']) # Spawning biomass becomes zooplankton
        else:
            dGonaddt = 0.
            dSomadt = 0.   
    
    
        #Feedback to NPZD2 model
        # Faeces and Pseudofaeces
        Fae = F[t] * ((par['epsilonP']*(1-par['AE_P'])*Phy[t])+ \
                     (par['epsilonZ']*(1-par['AE_Z'])*Zoo[t])+ \
                     (par['epsilonD']*(1-par['AE_D'])*SDet[t]))
        dLDetdt = dLDetdt + (Fae*(n_muss[t]/par['V']))     
                  
        # Remove eaten Phy, Zoo and SDet from water-column
        dPhydt =  dPhydt-((F[t] *par['epsilonP']*Phy[t])*(n_muss[t]/par['V']))
        dZoodt =  dZoodt-((F[t] *par['epsilonZ']*Zoo[t])*(n_muss[t]/par['V']))
        dSDetdt = dSDetdt -((F[t] *par['epsilonD']*SDet[t])*(n_muss[t]/par['V']))
        
        # Excretion into Ammonia
        dNH4dt = dNH4dt + (R[t]*n_muss[t]/par['V'])   
    
        # Population dynamics of mussels
        dn_mussdt = 0
    
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
        
        if Gonad[t+1] < 0:# If negative Gonad, take the biomass from Soma instead
            Soma[t+1] = Soma[t+1] + Gonad[t+1]
            Gonad[t+1] = 0. 
        
        B[t+1] = Soma[t+1] + Gonad[t+1]
        B_total[t+1] =  B[t+1] * n_muss[t+1] / par['V']
        n_muss[t+1] = n_muss[t] + (dn_mussdt * dt)    
            
        # Estimate Total Nitrogen
        TotN[t+1] = Phy[t+1] + Zoo[t+1] + SDet[t+1] + LDet[t+1] + NH4[t+1] + NO3[t+1] + ((B[t+1]*n_muss[t+1])/par['V'])
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
    output['Oxy'] = Oxy
    output['L_Temp'] = L_Temp 
    output['L_Salt'] = L_Salt 
    output['L_Oxy'] = L_Oxy
    output['L_Food'] = L_Food 
    output['B_total'] = B_total
    output['n_muss'] = n_muss


    print "Model run: DONE!!!"
    return output


    
def plot(output):
    '''
    Script to make plots
    '''
    # Import libraries
    import matplotlib.pyplot as plt
    
    # Plotting
    # Plankton NPZD2 ecosystem
    fig1, ax1 = plt.subplots(1,1,figsize=(15, 4))
    ax1.plot(output['time']/365,output['f_I'],'r-')
    ax1.plot(output['time']/365,output['mu'],'g-')
    ax1.plot(output['time']/365,output['L_NO3'],'b-')
    ax1.plot(output['time']/365,output['L_NH4'],'k-')
    ax1.set_xlabel('Time (years)')
    ax1.set_ylabel('Limitations (dimensionless)')
    ax1.set_title('PANEL 1: Fennel et al 2006 Model - Limitations')
    ax1.legend(['f_I','Mu','L_NO3','L_NH4'])
    plt.show()
    
    fig2, ax2 = plt.subplots(1,1,figsize=(15, 4))
    ax2.plot(output['time']/365,output['Phy'],'g-')
    ax2.plot(output['time']/365,output['Zoo'],'r-')
    ax2.plot(output['time']/365,output['SDet'],'k-')
    ax2.plot(output['time']/365,output['LDet'],'k-.')
    ax2.plot(output['time']/365,output['NH4'],'m-')
    ax2.plot(output['time']/365,output['NO3'],'c-')
    ax2.set_xlabel('Time (years)')
    ax2.set_ylabel('Nitrogen (mmol N m$^{-3}$)')
    ax2.set_title('PANEL 2: Fennel et al 2006 Model - State')
    plt.legend(['Phy','Zoo','SDet','LDet','NH4','NO3'])
    plt.show()
    
    fig3, ax3 = plt.subplots(1,1,figsize=(15, 4))
    ax3.plot(output['time']/365,output['Oxy'],'b-')
    ax3.set_xlabel('Time (years)')
    ax3.set_ylabel('Oxygen (mmol O$_2$ m$^{-3}$)')
    ax3.set_title('PANEL 3: Fennel et al 2006 Model - State - Oxygen')
    plt.legend(['Oxy'])
    plt.show()
    
    # Mussels - SHELL-E ecosystem
    fig4, ax4 = plt.subplots(1,1,figsize=(15, 4))
    ax4.plot(output['time']/365,output['L_Temp'],'r-') 
    ax4.plot(output['time']/365,output['L_Salt'],'c-') 
    ax4.plot(output['time']/365,output['L_Oxy'],'b')
    ax4.plot(output['time']/365,output['L_Food'],'g')
    ax4.set_ylabel('Limitations (dimensionless)')
    ax4.set_xlabel('Time (years)') 
    ax4.set_title('PANEL 4: SHELL-E Model - Limitations')
    ax4.legend(['L_Temp', 'L_Salt', 'L_Oxy','L_Food'])
    plt.show()

    fig5, ax5 = plt.subplots(1,1,figsize=(15, 4))
    ax5.plot(output['time']/365,output['B'],'r-') 
    ax5.plot(output['time']/365,output['Soma'],'b-') 
    ax5.plot(output['time']/365,output['Gonad'],'g.')
    ax5.set_ylabel('Nitrogen (mmol N m$^{-3}$)')
    ax5.set_xlabel('Time (years)')
    ax5.set_title('PANEL 5: SHELL-E Model - State')
    ax5.legend(['B', 'Soma', 'Gonad'])
    plt.show()
    
    # Total Nitrogen
    fig6, ax6 = plt.subplots(1,1,figsize=(15, 4))
    ax6.plot(output['time']/365,output['TotN'],'y.') 
    ax6.set_ylabel('Nitrogen (mmol N m$^{-3}$)')
    ax6.set_xlabel('Time (years)')
    ax6.set_title('PANEL 6: Total Nitrogen')
    ax6.legend(['TotN'])
    plt.show()
    return
    
if __name__ == '__main__':
    days, dt, par, InitCond = load_defaults()
    output = run(days,dt,InitCond,par)
    plot(output)