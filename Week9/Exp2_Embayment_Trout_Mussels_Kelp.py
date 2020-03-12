import model_NPZD2_SHELLE_openBoundary_withForcing_fishMacroalgae as model
import forcing
import plot_experiment as pexp
import pickle

experiment_name = 'Exp2_Embayment_fish'

days, dt, par, InitCond = model.load_defaults()

days = 365 * 1.5
forc = forcing.get(dt,days) # Load forcing

#------ Change the defaults ---------------
par['chi']  = 0.001 # No exchange of water in tank
par['X'] = 2000. # Basin length
par['Y'] = 100. # Basin width
par['Z'] = 10. # Basin depth
par['V'] = par['X'] * par['Y'] * par['Z']
InitCond['n_muss'] = 1e4
InitCond['n_Algae']= 1e6
# ----------------------------------------

Fish_levels = [1e3,5e3,7e3,1e4,2e4,5e4,1e5] # Units: number of mussels in tank (or bay)
#Fish_levels = [1e3,5e3,7e3,1e4,5e4,1e5] # Units: number of mussels in tank (or bay)
#Fish_levels = [1e3,5e3,7e3,1e4,5e4,7e4,1e5] # Units: number of mussels in tank (or bay)

multioutput = {}
for level in Fish_levels:
    InitCond['n_fish'] = float(level)
    multioutput[str(level)] = model.run(days,dt,InitCond,par,forc)
    
pickle.dump( multioutput, open( experiment_name + '.p', 'wb' ) )

print('Experiment is DONE!')

#%% Plotting ------------------------------------------
pexp.densityVSproductionFISH(multioutput,optimum=0.)
pexp.densityVSoxygenFISH(multioutput,optimum=0.)
pexp.densityVSammoniaFISH(multioutput,optimum=0.)