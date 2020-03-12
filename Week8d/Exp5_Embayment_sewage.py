import model_NPZD2_SHELLE_openBoundary_withForcing as model
import forcing
import plot_experiment as pexp
import pickle

experiment_name = 'Exp1_Tank_noFlow'

days, dt, par, InitCond = model.load_defaults()

days = 365 * 1.5
forc = forcing.get(dt,days) # Load forcing

#------ Change the defaults ---------------
par['chi']  = 0.001 # No exchange of water in tank
par['X'] = 2000. # Basin length
par['Y'] = 100. # Basin width
par['Z'] = 10. # Basin depth
par['V'] = par['X'] * par['Y'] * par['Z'] # <<<<<<<<<<<THIS IS NEW
par['sewage'] = 0.000001
# ----------------------------------------

#Mussel_levels = [1e4,1e5,1e6,1e7,1e8] # Units: number of mussels in tank (or bay)
#Mussel_levels = [2e9,5e9] # Units: number of mussels in tank (or bay)
Mussel_levels = [1e6,1e7,1e8,1e9,2e9,3e9,5e9,7e9,1e10] # Units: number of mussels in tank (or bay)

multioutput = {}
for level in Mussel_levels:
    InitCond['n_muss'] = float(level)
    multioutput[str(level)] = model.run(days,dt,InitCond,par,forc)
    
pickle.dump( multioutput, open( experiment_name + '.p', 'wb' ) )

print 'Experiment is DONE!'

#%% Plotting ------------------------------------------
pexp.densityVSproduction(multioutput,optimum=0.)
pexp.densityVSoxygen(multioutput,optimum=0.)
pexp.densityVSammonia(multioutput,optimum=0.)