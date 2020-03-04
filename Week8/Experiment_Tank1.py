import model_NPZD2_SHELLE_openBoundary_withForcing as model
import forcing
import plot_experiment as pexp
import pickle

days, dt, par, InitCond = model.load_defaults()

days = 365 * 1
forc = forcing.get(dt,days) # Load forcing

#------ Change the defaults ---------------
par['chi']  = 0.0 # No exchange of water in tank
par['X'] = 1. # Basin length
par['Y'] = 1. # Basin width
par['Z'] = 1. # Basin depth
par['V'] = par['X'] * par['Y'] * par['Z']
# ----------------------------------------

Mussel_levels = [0,1,2,4,8,16] # Units: number of mussels in tank (or bay)

Tank1_multioutput = {}
for level in Mussel_levels:
    InitCond['n_muss'] = float(level)
    Tank1_multioutput[str(level)] = model.run(days,dt,InitCond,par,forc)
    
pickle.dump( Tank1_multioutput, open( 'Tank1_multioutput.p', 'wb' ) )

print 'Experiment is DONE!'