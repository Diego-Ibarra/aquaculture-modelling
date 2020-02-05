'''
NPZD Adapted from:
Fennel et al (2006) Nitrogen cycling in the Middle Atlantic Bight: Results from a 
three-dimensional model and implications for the North Atlantic nitrogen budget.
GLOBAL BIOGEOCHEMICAL CYCLES, VOL. 20, GB3007, doi:10.1029/2005GB002456

Tilapia adapted from:
Yang Yi (1998) A bioenergetics growth model for Nile tilapia 
(Oreochromis niloticus) based on limiting nutrients and
fish standing crop in fertilized ponds. Aquaculture Engineering 18.157-173

Oysters adapted from:
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
    days =365 * 3 # Three year
    dt   = 0.01 # units: days    
    
    # Parameters you MUST change to reflect your site & conditions ========================
    par = {}    
    
    # Physical characteristics of embayment
    par['chi'] = 0.00 # Exchange rate between pond and open ocean (d-1)
    par['X'] = 2000.#2000 # Basin length
    par['Y'] = 100.#200 # Basin width
    par['Z'] = 10. # Basin depth
    par['V'] = par['X'] * par['Y'] * par['Z']
    
    #Parameters Oysters
    par['lamda_nat'] = 0.0 # 0.00137 (d-1)
    par['lamda_harvest'] = 0.0 # 0.001 (d-1)
    # Parameters Tilapia
    par['feeding_rate'] = 0.01 # Fraction of body weight feed per day (d-1)
    par['lamda_nat_tilapia'] = 0  # Mortality rate tilapia (d-1)

    
    # Initial conditions you MUST change to reflect your site & conditions ========================
    InitCond = {}
    # Oyster parameters
    InitCond['conc_oys'] = 2.
    # Tilapia
    InitCond['W'] = 13. # (g)
    InitCond['n_tilapia'] = 10000.  # Number of tilapia in pond (i.e inmortal tilapia)
    
    
    
    
    
    
    
    
    
    
    # DON'T CHANGE THESE !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    # Other Parameters - DON'T CHANGE THESE !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    # Plankton parameters
    par['mu0']   = 0.69  
    par['kNO3']  = 0.5    
    par['kNH4']  = 0.5  
    par['alpha'] = 0.125  
    par['gmax']  = 0.6  #Original 0.6
    par['kP']    = 2
    par['mP']    = 0.15    
    par['tau']   = 0.005 
    par['thetaMax']= 0.053
    par['beta']  = 0.75 
    par['lBM']   = 0.01   
    par['lE']    = 0.01
    par['mZ']    = 0.025
    par['rSD']   = 0.3 # Original 0.03
    par['rLD']   = 0.1 # # Original 0.01
    par['nmax']  = 0.05
    par['kI']    = 0.1
    par['I0']    = 0.0095
    par['wP']    = 0.1
    par['wS']    = 0.1
    par['wL']    = 1.
    # Oyster Parameters
    par['AE_P']  = 0.9  
    par['AE_D']  = 0.2    
    par['AE_Z']  = 0.3  
    par['Bpub']  = 0.43  
    par['Fmax_ref']= 0.03#0.025
    par['GT']    = 0.44
    par['KTempH']= 0.1    
    par['KTempL']= 0.5 
    par['KSaltL']= 0.25
    par['KOxyL'] = 0.02 
    par['KFood'] = 1.    
    par['KRE']   = 0.86
    par['OxyL']  = 17.5
    par['Rm']    = 0.002
    par['SaltL'] = 10.
    par['TempH'] = 35.
    par['TempL'] = 10.
    par['beta']  = 0.12
    par['epsilonP'] = 1.
    par['epsilonD'] = 0.5
    par['epsilonZ'] = 0.3

    #Tilapia
    par['b'] = 0.62 # efficiency of food assimilation(dimensionless)
    par['a'] = 0.53 # fraction of food assimilated used for feeding catabolism(dimensionless)
    par['h'] = 0.8 #(dimensionless)
    par['s'] = 1#17.31 #(dimensionless)
    par['m'] = 0.67 #(dimensionless)
    par['n'] = 0.81 #(dimensionless)
    par['PNPP/B'] = 1 # Potential net primary production (PNPP) to standing crop (B) of Tilapia (dimensionless)
    par['kmin'] = 0.00133
    par['j'] = 0.0132
    par['T'] = 35 # Temperature (Degrees C)
    par['Topt'] = 33 # Minimum Temperature for survival (Degrees C)
    par['Tmin'] = 15 # Minimum Temperature for survival (Degrees C)
    par['Tmax'] = 41 # Minimum Temperature for survival (Degrees C)
    par['DOmin'] = 9.375 #Minimum Dissolved Oxygen - mmol O2 / M3
    par['DOcrit'] = 32.25 # Critical Dissolved Oxygen - mmol O2 / M3
    par['UIA'] = 0.07 # Unionized Ammonia Concentration mg/l
    par['UIAmax'] = 77.8  # maximum Unionized Ammonia Concentration - mmol N m3
    par['UIAcrit'] = 3.333 # critical Unionized Ammonia Concentration mmol N m3
    par['WW2Nitrogen'] = 3.085 #Convertion factor 3.085 mmol N / g Fish
    
    # Other physical parameters (for Oxygen exchange with Atmosphere)
    par['uwind'] = 0.5
    par['vwind'] = 0.5
    
    # Other Initial Conditions
    InitCond['Phy'] = 0.5
    InitCond['Zoo'] = 0.2
    InitCond['SDet'] = 0.5 
    InitCond['LDet'] = 0.3
    InitCond['NH4'] = 0.1
    InitCond['NO3'] = 5.
    InitCond['Oxy'] = 340. #Oxygen 
    InitCond['Soma'] = 1.0
    InitCond['Gonad'] = 0.0
    InitCond['n_oys'] = InitCond['conc_oys'] * par['V']
    
    return days, dt, par, InitCond
    



    
def run_model(days,dt,InitCond,par,forc):
    '''
    This is your model. Do a brief description.

    '''
    # Import libraries
    import numpy as np
    
    print 'Starting model run with ' + str(InitCond['conc_oys']) + ' oysters/m3 ...'
    print '...and ' + str(InitCond['n_tilapia'] * InitCond['W'] / par['V']) + ' tilapia/m3 ...'
    print '...fed artificial food at a rate of ' + str(par['feeding_rate']) + ' percent-body-weight per day'
    
    # Make sure n_oys is correctly estimated
    InitCond['n_oys'] = InitCond['conc_oys'] * par['V']
    
    # Setup the framework (calculate timestemps, create zero vectors, create time vector)
    NoSTEPS = int(days / dt) # Calculates the number of steps by dividing days by dt and rounding down
    time = np.linspace(0,days,NoSTEPS) # Makes and vector array of equally spaced numbers from zero to "days"
    
    # Create zero-vectors
    Phy = np.zeros((NoSTEPS,),float) # makes a vector array of zeros (size: NoSTEPS rows by ONE column)
    Zoo = np.zeros((NoSTEPS,),float) # same as above
    SDet = np.zeros((NoSTEPS,),float) # Biomass - same as above 
    LDet = np.zeros((NoSTEPS,),float) # same as above
    NH4 = np.zeros((NoSTEPS,),float) # same as above
    NO3 = np.zeros((NoSTEPS,),float) # same as above
    Oxy = np.zeros((NoSTEPS,),float) # same as above

    
    mu = np.zeros((NoSTEPS,),float) # same as above
    f_I = np.zeros((NoSTEPS,),float) # same as above
    L_NO3 = np.zeros((NoSTEPS,),float) # same as above
    L_NH4 = np.zeros((NoSTEPS,),float) # same as above
    airwater_O2_flux = np.zeros((NoSTEPS,),float) # same as above
    TotN = np.zeros((NoSTEPS,),float) # same as above
    
    # Oysters
    Soma = np.zeros((NoSTEPS,),float) # makes a vector array of zeros (size: NoSTEPS rows by ONE column)
    Gonad = np.zeros((NoSTEPS,),float) # same as above
    B = np.zeros((NoSTEPS,),float) # Biomass - same as above 
    B_conc = np.zeros((NoSTEPS,),float) # Biomass - same as above 
    L_Temp = np.zeros((NoSTEPS,),float) # same as above
    L_Salt = np.zeros((NoSTEPS,),float) # same as above
    L_Oxy = np.zeros((NoSTEPS,),float) # same as above
    L_Food = np.zeros((NoSTEPS,),float) # same as above
    F = np.zeros((NoSTEPS,),float) # same as above
    A = np.zeros((NoSTEPS,),float) # same as above
    R = np.zeros((NoSTEPS,),float) # same as above
    RE = np.zeros((NoSTEPS,),float) # same as above
    Spawning = np.zeros((NoSTEPS,),float) # same as above
    n_oys = np.zeros((NoSTEPS,),float) # same as above
    CumulativeHarvest = np.zeros((NoSTEPS,),float) # same as above
    
    # Tilapia
    W = np.zeros((NoSTEPS,),float)
    n_tilapia= np.zeros((NoSTEPS,),float)
    delta= np.zeros((NoSTEPS,),float)
    upsilon= np.zeros((NoSTEPS,),float)
    tauT= np.zeros((NoSTEPS,),float)
    f= np.zeros((NoSTEPS,),float)

    
    
    # Initializing with initial conditions
    #NPZD
    Phy[0] = InitCond['Phy']
    Zoo[0] = InitCond['Zoo']
    SDet[0] = InitCond['SDet']
    LDet[0] = InitCond['LDet']
    NH4[0] = InitCond['NH4']
    NO3[0] = InitCond['NO3']
    Oxy[0] = InitCond['Oxy']
    # Oysters
    Soma[0] = InitCond['Soma']
    Gonad[0] = InitCond['Soma']
    B[0] = InitCond['Soma'] + InitCond['Gonad']
    Spawning[0] = 0.
    n_oys[0] = InitCond['n_oys']
    B_conc[0] = B[0] * n_oys[0] / par['V']
    # Tilapia
    W[0] = InitCond['W']
    n_tilapia[0] = InitCond['n_tilapia']
    
    # Forcing
    Salt = forc['Salt']
    Temp = forc['Temp']
    I = forc['I']
    Phy0 = forc['Phy']
    Zoo0 = forc['Zoo']
    SDet0 = forc['SDet']
    LDet0 = forc['LDet']
    NH40 = forc['NH4']
    NO30 = forc['NO3']
    Oxy0 = forc['Oxy']


    # *****************************************************************************
    # MAIN MODEL LOOP *************************************************************
    for t in range(0,NoSTEPS-1):
        muMax = par['mu0'] * (1.066 ** Temp[t]) # text
        
        f_I[t] = (par['alpha']*I[t])/(np.sqrt(muMax**2+((par['alpha']**2)*(I[t]**2)))) #Eq5
        
        L_NO3[t] = (NO3[t]/(par['kNO3']+NO3[t])) * (1/(1+(NH4[t]/par['kNH4']))) #Eq3
        
        L_NH4[t] = NH4[t]/(par['kNH4']+NH4[t]) # Eq 4
    
        mu[t] =muMax * f_I[t] * (L_NO3[t] + L_NH4[t]) # Eq2
    
        g = par['gmax'] * ((Phy[t]**2)/(par['kP']+(Phy[t]**2)))
        
        n = par['nmax'] * (1 - max(0,(I[t]-par['I0'])/(par['kI']+I[t]-par['I0'])))
        
        n_O2 = (Oxy[t]/(3.+Oxy[t]))
    
        dPhydt = (mu[t] * Phy[t]) - \
                 (g  * Zoo[t]) - \
                 (par['mP'] * Phy[t]) - \
                 (par['tau']*(SDet[t]+Phy[t])*Phy[t]) - \
                 (par['wP']*Phy[t]/par['Z']) # Eq1
                 
        dZoodt = (g * par['beta'] * Zoo[t]) - \
                 (par['lBM']*Zoo[t]) - \
                 (par['lE']*((Phy[t]**2)/(par['kP']+(Phy[t]**2)))*par['beta']*Zoo[t]) - \
                 (par['mZ']*(Zoo[t]**2))#Eq10
        
        dSDetdt = (g * (1-par['beta']) * Zoo[t]) + \
                  (par['mZ']*(Zoo[t]**2)) + \
                  (par['mP'] * Phy[t]) - \
                  (par['tau']*(SDet[t]+Phy[t])*SDet[t]) - \
                  (par['rSD']*SDet[t]) - \
                  (par['wS']*SDet[t]/par['Z'])
                  
        dLDetdt = (par['tau']*((SDet[t]+Phy[t])**2)) - \
                  (par['rLD']*LDet[t]) - \
                  (par['wL']*LDet[t]/par['Z'])
                  
        dNO3dt = -(muMax * f_I[t] * L_NO3[t] * Phy[t]) + \
                  (n * n_O2 * NH4[t])
                 
        dNH4dt = -(muMax * f_I[t] * L_NH4[t] * Phy[t]) - \
                  (n * n_O2 * NH4[t]) + \
                  (par['lBM'] * Zoo[t]) + \
                  (par['lE']*((Phy[t]**2)/(par['kP']+(Phy[t]**2)))*par['beta']*Zoo[t]) + \
                  (par['rSD']*SDet[t]) + \
                  (par['rLD']*LDet[t]) + \
                  (par['wP']*Phy[t]/par['Z']) + \
                  (par['wS']*SDet[t]/par['Z']) + \
                  (par['wL']*LDet[t]/par['Z'])


        # OYSTERS -------------------------------------------------------------
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
        if n_oys[t] == 0.: 
            Spawning[t] = 0.
            dGonaddt = 0
            dSomadt = 0.
        elif Gonad[t]/B[t] < par['GT']:
            Spawning[t] = 0.
            dGonaddt = (A[t]-R[t]) * RE[t]
            dSomadt =  (A[t]-R[t]) * (1.-RE[t])
            offset = Gonad[t] + (dGonaddt * dt)
            if offset < 0 : # If Gonad is going to be negative... don't apply dynamic allocation
                dGonaddt = 0.
                dSomadt = A[t]-R[t]
        elif Gonad[t]/B[t] >= par['GT']:         
            Spawning[t] = Gonad[t]
            dGonaddt = 0.
            dSomadt = A[t]-R[t]

        #Feedback to NPZD2 model
        # Faeces and Pseudofaeces
        Fae = F[t] * ((par['epsilonP']*(1-par['AE_P'])*Phy[t])+ \
                     (par['epsilonZ']*(1-par['AE_Z'])*Zoo[t])+ \
                     (par['epsilonD']*(1-par['AE_D'])*SDet[t]))
                     
        dLDetdt = dLDetdt + (Fae*(n_oys[t]/par['V']))    
                  
        # Remove eaten Phy, Zoo and SDet from water-column
        dPhydt =  dPhydt-((F[t] *par['epsilonP']*Phy[t])*(n_oys[t]/par['V']))
        dZoodt =  dZoodt-((F[t] *par['epsilonZ']*Zoo[t])*(n_oys[t]/par['V']))
        dSDetdt = dSDetdt -((F[t] *par['epsilonD']*SDet[t])*(n_oys[t]/par['V']))

        # Oysters population
        Harvest = par['lamda_harvest'] * n_oys[t]
        NatMortality = par['lamda_nat'] * n_oys[t]
        dn_oysdt = -NatMortality - Harvest
   
        # Excretion into Ammonia
        dNH4dt = dNH4dt + ((R[t]*n_oys[t]/par['V']) + ((NatMortality*B[t])/par['V']))       
      
        
        
        # Tilapia Sub-model ---------------------------------------------------
        if Oxy[t] > par['DOcrit']:
            delta[t] = 1
        elif par['DOmin'] <= Oxy[t] and Oxy[t] <= par['DOcrit']:
            delta[t] = (Oxy[t]-par['DOmin'])/(par['DOcrit']-par['DOmin'])
        elif Oxy[t] < par['DOmin']:
            delta[t] = 0
        # Equation 3  (a,b,c) respectively
        if NH4[t] < par['UIAcrit']:      
            upsilon[t] = 1 # !!!TO DO
        elif par['UIAcrit'] <= NH4[t]  and NH4[t]  <= par['UIAmax']:
            upsilon[t] = (par['UIAmax']-NH4[t] )/(par['UIAmax'] - par['UIAcrit'])
        elif NH4[t]  > par['UIAmax']:
            upsilon[t] = 0
               
        if Temp[t] < par['Topt']:
            tauT[t] = np.exp(-4.6 * (((par['Topt']-Temp[t])/(par['Topt'] - par['Tmin']))**4))
        elif Temp[t] >= par['Topt']:
            tauT[t] = np.exp(-4.6 * (((Temp[t]-par['Topt'])/(par['Tmax'] - par['Topt']))**4))
            
        # Equation 4 (a,b,c) respectively
        k = par['kmin'] * np.exp(par['j']*(Temp[t] - par['Tmin']))  # Equation 12
        
        f[t] = 1 - np.exp(-par['s']*(Phy[t]/(W[t]*n_tilapia[t]/par['V']))) # Equation 9
        
        dRdt_Phy =  tauT[t] * delta[t] * upsilon[t] * par['h'] * f[t] * (W[t]**par['m'])# Phytoplanklton Consumption rate
        
        dRdt_Food =  tauT[t] * delta[t] * upsilon[t] * par['feeding_rate'] * W[t]# Food Consumption rate
        
        dRdt = dRdt_Phy + dRdt_Food
        
        FoodAbsorbed = par['b'] * dRdt #Deducted From Eq2
        
        Faeces = (1 - par['b']) * dRdt # Deducted
        
        Resp = (par['a'] * FoodAbsorbed) + (k*(W[t]**par['n'])) #Deducted From Eq2
                
        dWdt = FoodAbsorbed - Resp#Equation 2 - Final Growth Rate

        # Oysters population
        NatMortality_Tilapia = par['lamda_nat_tilapia'] * n_tilapia[t]  
        dn_tilapiadt = -NatMortality_Tilapia
        
        #Feedback to NPZD2 model
        # Faeces
        dLDetdt = dLDetdt + (Faeces * n_tilapia[t] / par['V'])
                  
        # Remove eaten Phy, Zoo and SDet from water-column
        dPhydt =  dPhydt - (dRdt_Phy * n_tilapia[t] / par['V'])
        
        # Excretion into Ammonia
        dNH4dt = dNH4dt + (Resp * n_tilapia[t] / par['V']) + ((NatMortality_Tilapia*W[t])/par['V'])
        # End of Tilapia ------------------------------------------------------        
        
        
        
        
        
        
        # Oxygen sub-model =========================================
        
        # Parameters
        OA0 = 2.00907       # Oxygen
        OA1 = 3.22014       # saturation
        OA2 = 4.05010       # coefficients
        OA3 = 4.94457
        OA4 =-0.256847
        OA5 = 3.88767
        OB0 =-0.00624523
        OB1 =-0.00737614
        OB2 =-0.0103410
        OB3 =-0.00817083
        OC0 =-0.000000488682
        rOxNO3= 8.625       # 138/16
        rOxNH4= 6.625       # 106/16
        l2mol = 1000.0/22.9316 # liter to mol
        
        #-----------------------------------------------------------------------
        #  Surface O2 gas exchange.
        #-----------------------------------------------------------------------
        
        #  Compute surface O2 gas exchange.
        cff2=0.31*(24.0/100.0)
        
        #  Compute O2 transfer velocity : u10squared (u10 in m/s)
        u10squ=(par['uwind']*par['uwind'])+(par['vwind']*par['vwind'])
        
        # Calculate the Schmidt number for O2 in sea water (Wanninkhof, 1992).
        SchmidtN_Ox=1953.4-Temp[t]*(128.0-Temp[t]*(3.9918-Temp[t]*0.050091))
        cff3=cff2*u10squ*np.sqrt(660.0/SchmidtN_Ox)        
        
        #  Calculate O2 saturation concentration using Garcia and Gordon
        #  L&O (1992) formula, (EXP(AA) is in ml/l).        
        TS=np.log((298.15-Temp[t])/(273.15+Temp[t]))        
        
        AA=OA0+TS*(OA1+TS*(OA2+TS*(OA3+TS*(OA4+TS*OA5))))+ \
           Salt[t]*(OB0+TS*(OB1+TS*(OB2+TS*OB3)))+ \
           OC0*Salt[t]*Salt[t]
        
        # Convert from ml/l to mmol/m3.
        O2satu=l2mol*np.exp(AA)        
        
        # Add in O2 gas exchange.
        O2_Flux = cff3*(O2satu-Oxy[t])
        
        airwater_O2_flux[t] = O2_Flux * (1./par['Z'])
        
        dOxydt = airwater_O2_flux[t]
        
        
        # Production via Photosynthesys
        dOxydt = dOxydt + (muMax * f_I[t] * L_NO3[t] * Phy[t] * rOxNO3) # New production
        dOxydt = dOxydt + (muMax * f_I[t] * L_NH4[t] * Phy[t] * rOxNH4) # Regenerated production
        
        # Respiration
        dOxydt = dOxydt - (((par['lBM']*Zoo[t]) - \
                           (par['lE']*((Phy[t]**2)/(par['kP']+(Phy[t]**2)))*par['beta']*Zoo[t]) - \
                           (par['mZ']*(Zoo[t]**2))) * rOxNH4) # Zooplankton
 
        dOxydt = dOxydt - ((n * n_O2 * NH4[t])* rOxNH4 * 2) # Nitrification 
 
        dOxydt = dOxydt - (((par['rSD']*SDet[t])+(par['rLD']*LDet[t])) * rOxNH4) #S and L Detritus remineralization
        
        dOxydt = dOxydt - (((par['wS']*SDet[t]/par['Z']) + \
                          (par['wL']*LDet[t]/par['Z'])) * rOxNH4) #S and L Detritus remineralization in sediments
        
        dOxydt = dOxydt - ((((R[t]*n_oys[t])/par['V']) + ((NatMortality*B[t])/par['V'])) * rOxNH4)  #oysters
        dOxydt = dOxydt - (((Resp * n_tilapia[t] / par['V']) + ((NatMortality_Tilapia*W[t])/par['V'])) * rOxNH4)  #Tilapia Respiration
       


        
        # Physical Model ======================================================
        dPhydt = dPhydt + (par['chi'] * (Phy0[t] - Phy[t]))
        dZoodt = dZoodt + (par['chi'] * (Zoo0[t] - Zoo[t]))
        dNH4dt = dNH4dt + (par['chi'] * (NH40[t] - NH4[t]))
        dNO3dt = dNO3dt + (par['chi'] * (NO30[t] - NO3[t]))
        dSDetdt = dSDetdt + (par['chi'] * (SDet0[t] - SDet[t]))
        dLDetdt = dLDetdt + (par['chi'] * (LDet0[t] - LDet[t]))
        dOxydt = dOxydt + (par['chi'] * (Oxy0[t] - Oxy[t]))
                
        

        # Update and step ----------------------------------------------------
        Phy[t+1]  = Phy[t]  + (dPhydt * dt)
        Zoo[t+1]  = Zoo[t]  + (dZoodt * dt) + ((Spawning[t]*n_oys[t])/par['V'])
        SDet[t+1] = SDet[t] + (dSDetdt * dt)
        LDet[t+1] = LDet[t] + (dLDetdt * dt)
        Oxy[t+1]  = max(0,Oxy[t] +  (dOxydt * dt))
        NH4[t+1]  = NH4[t]  + (dNH4dt * dt)
        NO3[t+1]  = NO3[t] +  (dNO3dt * dt)
        if NO3[t+1] <= 0.001:
            offset = NO3[t+1]
            NH4[t+1] = NH4[t+1] + offset
            NO3[t+1] = NO3[t+1] - offset
        # Oysters
        Soma[t+1] = Soma[t] + (dSomadt * dt)
        Gonad[t+1] = Gonad[t] + (dGonaddt * dt) - Spawning[t]
        B[t+1] = Soma[t+1] + Gonad[t+1]
        n_oys[t+1] = max(0,n_oys[t] + (dn_oysdt * dt))
        CumulativeHarvest[t+1] = CumulativeHarvest[t] + Harvest
        B_conc[t+1] =  B[t+1] * n_oys[t+1] / par['V']
        # Tilapia
        W[t+1] = W[t] + (dWdt*dt)
        n_tilapia[t+1] = max(0,n_tilapia[t] + (dn_tilapiadt * dt))
        # Estimate Total Nitrogen
        TotN[t+1] = Phy[t+1] + Zoo[t+1] + SDet[t+1] + LDet[t+1] + NH4[t+1] + \
                    NO3[t+1] + ((B[t+1]*n_oys[t+1])/par['V']) + ((W[t+1]*n_tilapia[t+1])/par['V'])
    # end of main model LOOP*******************************************************
    # *****************************************************************************

    # Pack output into dictionary
    output = {}
    output['par'] = par
    output['InitCond'] = InitCond
    output['time'] = time
    output['Phy'] = Phy
    output['Zoo'] = Zoo
    output['SDet'] = SDet
    output['LDet'] = LDet
    output['NH4'] = NH4
    output['NO3'] = NO3
    output['Oxy'] = Oxy
    output['I'] = I
    output['mu'] = mu
    output['f_I'] = f_I
    output['L_NO3'] = L_NO3
    output['L_NH4'] = L_NH4
    output['TotN'] = TotN
    output['airwater_O2_flux'] = airwater_O2_flux
    output['Soma'] = Soma
    output['Gonad'] = Gonad
    output['B'] = B
    output['n_oys'] = n_oys
    output['Spawning'] = Spawning
    output['CumulativeHarvest'] = CumulativeHarvest
    output['F'] = F
    output['L_Temp'] = L_Temp
    output['L_Salt'] = L_Salt
    output['L_Oxy'] = L_Oxy
    output['L_Food'] = L_Food
    output['B_conc'] = B_conc
    #Tilapia
    output['W'] = W
    output['n_tilapia'] = n_tilapia
    output['delta'] = delta
    output['upsilon'] = upsilon
    output['tauT'] = tauT
    output['f'] = f

    print "Model run: DONE!!!"
    return output


    
def plot_model(output):
    '''
    Script to make plots
    '''
    # Import libraries
    import matplotlib.pyplot as plt
    
    # Plotting
    fig, (ax, ax2, ax3) = plt.subplots(3,1,figsize=(13,13))
    ax.plot(output['time']/365,output['Phy'],'g-')
    ax.plot(output['time']/365,output['Zoo'],'r-')
    ax.plot(output['time']/365,output['SDet'],'k-')
    ax.plot(output['time']/365,output['LDet'],'k-.')
    ax.plot(output['time']/365,output['NH4'],'m-')
    ax.plot(output['time']/365,output['NO3'],'c-')
    ax.plot(output['time']/365,output['B_conc'],'r.')
#    ax.plot(output['time']/365,output['W']*output['n_tilapia']*output['par']['V'],'r.')
    ax.set_ylabel('Nitrogen \n (mmol N m$^{-3}$)')
    ax.set_title('Ecosystem Model - Plankton Ecosystem')
    ax.legend(['Phy','Zoo','SDet','LDet','NH4','NO3','B_conc'])
    
    ax2.plot(output['time']/365,output['Oxy'],'b-')
    ax2.set_ylabel('Oxygen \n (mmol O2 m$^{-3}$)')
    ax2.legend(['f_I','Mu','L_NO3','L_NH4'])
    
    ax3.plot(output['time']/365,output['f_I'],'r-')
    ax3.plot(output['time']/365,output['mu'],'g-')
    ax3.plot(output['time']/365,output['L_NO3'],'b-')
    ax3.plot(output['time']/365,output['L_NH4'],'k-')
    ax3.set_ylabel('Plankton Diagnostics \n (dimensionless)')
    ax3.legend(['f_I','Mu','L_NO3','L_NH4'])


    
    fig2, (ax, ax2, ax3, ax4) = plt.subplots(4,1,figsize=(13,13))
    ax.plot(output['time']/365,output['B'],'r.')
    ax.plot(output['time']/365,output['Gonad'],'b-')
    ax.plot(output['time']/365,output['Soma'],'k-')
    ax.set_ylabel('Nitrogen \n (mmol N m$^{-3}$)')
    ax.legend(['B','Gonad','Soma'])
    ax.set_title('Ecosystem Model - oysters')
    
    ax2.plot(output['time']/365,output['n_oys'],'g-')
    ax2.set_ylabel('Number of oysters in bay')
    ax2.legend(['Total Number of \n oysters in bay'])

    ax3.plot(output['time']/365,output['F']*1000/24,'g-')
    ax3.set_ylabel('Filtration rate \n (L ind$^{-1}$ h$^{-1}$)')
    ax3.legend(['F'])

    ax4.plot(output['time']/365,output['L_Temp'],'b-')
    ax4.plot(output['time']/365,output['L_Salt'],'m-')
    ax4.plot(output['time']/365,output['L_Oxy'],'k-')
    ax4.plot(output['time']/365,output['L_Food'],'r-')
    ax4.legend(['L_Temp','L_Salt','L_Oxy','L_Food'])
    ax4.set_ylabel('Oysters Diagnostics \n (dimensionless)')
    ax4.set_xlabel('Time (years)')


    fig2, (ax, ax2, ax3, ax4) = plt.subplots(4,1,figsize=(13,13))
    ax.plot(output['time']/365,output['W'],'r.')
    ax.set_ylabel('Nitrogen \n (mmol N m$^{-3}$)')
    ax.legend(['Nitrogen Tilapia'])
    ax.set_title('Ecosystem Model - Tilapia')

    ax2.plot(output['time']/365,output['W']/output['par']['WW2Nitrogen'],'r.')
    ax2.set_ylabel('Wet weight \n (g)')
    ax2.legend(['Wet weight Tilapia'])

    ax3.plot(output['time']/365,output['n_tilapia'],'g-')
    ax3.set_ylabel('Number of Tilapia in bay')
    ax3.legend(['Total Number of \n tilapia in bay'])

    ax4.plot(output['time']/365,output['delta'],'b-')
    ax4.plot(output['time']/365,output['upsilon'],'m-')
    ax4.plot(output['time']/365,output['tauT'],'k-')
    ax4.plot(output['time']/365,output['f'],'r-')
    ax4.legend(['delta (O2)','upsilon (NH4)','TauT (Temp)','f (Phy)'])
    ax4.set_ylabel('Tilapia Diagnostics \n (dimensionless)')
    ax4.set_xlabel('Time (years)')

    fig3, (ax) = plt.subplots(1,1)
    ax.plot(output['time']/365,output['TotN'],'y.')
    ax.legend(['TotN'])
    ax.set_title('Ecosystem Model - Total Nitrogen')
    ax.set_ylabel('Nitrogen \n (mmol N)')
    ax.set_xlabel('Time (years)')

    plt.show()
    return

    
    
    
if __name__ == '__main__':
    import new_load_forcing

    days, dt, par, InitCond = load_defaults()
    forc = new_load_forcing.get_forcing(dt,days)
    output = run_model(days,dt,InitCond,par,forc)
#    load_forcing.plot_forcing(dt,days,forc)
    plot_model(output)
