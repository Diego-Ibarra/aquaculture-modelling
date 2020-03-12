import model_NPZD2_SHELLE_openBoundary_withForcing as model
import forcing
import plot_experiment as pexp
import pickle

days, dt, par, InitCond = model.load_defaults()

days = 365 * 1.5
forc = forcing.get(dt,days) # Load forcing

#------ Change the defaults ---------------
par['chi']  = 0.001 # No exchange of water in tank
par['X'] = 2000 # Basin length
par['Y'] = 100 # Basin width
par['Z'] = 10 # Basin depth
par['V'] = par['X'] * par['Y'] * par['Z']
# ----------------------------------------

Mussel_levels = [10000,100000,1000000,5000000,10000000,50000000,100000000,100000000] # Units: number of mussels in tank (or bay)

Tank1_multioutput = {}
for level in Mussel_levels:
    InitCond['n_muss'] = float(level)
    Tank1_multioutput[str(level)] = model.run(days,dt,InitCond,par,forc)
    
pickle.dump( Tank1_multioutput, open( 'Tank1_multioutput.p', 'wb' ) )

print 'Experiment is DONE!'


pexp.densityVSproduction(Tank1_multioutput)
#fig1.savefig('Exp4_Embayment_HigherFlow_densVsprod.png')

pexp.densityVSoxygen(Tank1_multioutput)
#fig2.savefig('Exp4_Embayment_HigherFlow_densVsoxy.png')

pexp.densityVSammonia(Tank1_multioutput)
#fig3.savefig('Exp4_Embayment_HigherFlow_densVsammonia.png')
