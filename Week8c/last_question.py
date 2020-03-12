import model_NPZD2_SHELLE_openBoundary_withForcing as model
import forcing
import plot_experiment as pexp

experiment_name = 'Exp1_Tank_noFlow'

days, dt, par, InitCond = model.load_defaults()

days = 365 * 1.5
forc = forcing.get(dt,days) # Load forcing

#------ Change the defaults ---------------
par['chi']  = 0.2 # No exchange of water in tank
par['X'] = 1. # Basin length
par['Y'] = 1. # Basin width
par['Z'] = 1. # Basin depth
par['V'] = par['X'] * par['Y'] * par['Z'] # <<<<<<<<<<<THIS IS NEW
# ----------------------------------------

Mussel_levels = [0,5,10,20,30,50,100,150] # Units: number of mussels in tank (or bay)

multioutput = {}
for level in Mussel_levels:
    InitCond['n_muss'] = float(level)
    multioutput[str(level)] = model.run(days,dt,InitCond,par,forc)

print 'Experiment is DONE!'

pexp.densityVSproduction(multioutput,optimum=0.)