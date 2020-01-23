def load_defaults():
    print('Loading defaults...')
    # Framework
    days = 365 * 3 # One year
    dt   = 0.01 # units: days    
    
    # Parameters
    par = {}
    
    # Initial conditions
    InitCond = {}
    return  days, dt, par, InitCond
    
def run(days, dt, par, InitCond):
    print('Running model...')
    # Import libraries
    import numpy as np
    
    # Setup the framework 
    NoSTEPS = int(days / dt) # Calculates the number of steps 
    time = np.linspace(0,days,NoSTEPS) # Makes vector array of equally spaced numbers 
    
    # Create arrays of zeros
    
    # Initializing with initial conditions
    
    # *****************************************************************************
    # MAIN MODEL LOOP *************************************************************
    for t in range(0,NoSTEPS-1):
        a = 0 #DUMMY LINE
        # Update and step ------------------------------
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