import model_NPZD2_SHELLE_embayment_O2 as model
reload(model)

days, dt, par, InitCond = model.load_defaults()

par['X'] = 2000.# Basin length
par['Y'] = 100. # Basin width
par['Z'] = 10.  # Basin depth
par['V'] = par['X'] * par['Y'] * par['Z']

InitCond['n_muss'] = 4 * 10**7. # number of mussels in basin    
output = model.run(days,dt,InitCond,par)


InitCond['n_muss'] = 4 * 10**4. # number of mussels in basin      
output2 = model.run(days,dt,InitCond,par)


#%% PLotting
import matplotlib.pyplot as plt
fig, (ax) = plt.subplots(1,1,figsize=(15, 4))
ax.plot(output['n_muss'][0],output['Oxy'].mean(),'b.')
ax.plot(output2['n_muss'][0],output2['Oxy'].mean(),'b.')


