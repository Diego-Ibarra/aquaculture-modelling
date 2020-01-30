def load_defaults():
    print('Loading defaults...')
    # Framework
    days = 365 * 3 # One year
    dt   = 0.01 # units: days    
    
    # Parameters
    par = {}
    par['mZ']=1
    # EXAMPLE: par['mu'] = 1
    
    # Initial conditions
    InitCond = {}
    # EXAMPLE: InitCond['P'] = 1
    
    return  days, dt, par, InitCond
    
def run(days, dt, par, InitCond):
    print('Running model...')
    # Import libraries
    import numpy as np
    
    # Setup the framework 
    NoSTEPS = int(days / dt) # Calculates the number of steps 
    time = np.linspace(0,days,NoSTEPS) # Makes vector array of equally spaced numbers 
    
    # Create arrays of zeros
    # EXAMPLE: P = np.zeros((NoSTEPS), float)
    Z = np.zeros((NoSTEPS), float)
    
    # Initializing with initial conditions
    # EXAMPLE: P[0] = InitCond['P']
    Z[0] = 1 #Todo later
    
    # *****************************************************************************
    # MAIN MODEL LOOP *************************************************************
    for t in range(0,NoSTEPS-1):
        # Estimate limiting functions
        
        #Estimate processes
        a = 0 #DUMMY LINE
        
        #Estimate Time rate of change of all State Variables (dXdt) ----
        #EXAMPLE: dPdt = 0
        dZdt = -par['mZ']* Z[t]
        
        # Update and step (time-stepping) ------------------------------
        # EXAMPLE: P[t+1] = P[t] + (dPdt * dt)
        Z[t+1] = Z[t] + (dZdt * dt)
        
    # end of main model LOOP*******************************************************
    # *****************************************************************************

    # Pack output into dictionary
    output = {}
    output['time'] = time

    print('Model run: DONE!!!')
    return  output

def plot(output):
    import matplotlib.pyplot as plt 
    # Plotting                      
    fig, (ax) = plt.subplots(1,1)   
    plt.show()                      
    return

if __name__ == "__main__":
    print('Executing my_module.py')
    print('--------------------')
    
    days, dt, par, InitCond = load_defaults()
    output = run(days, dt, par, InitCond)
    plot(output)
    
    print('--------------------')