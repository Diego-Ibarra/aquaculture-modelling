# -*- coding: utf-8 -*-
"""
Created on Tue Oct 03 20:45:01 2017

@author: cerc-user
"""

#%% Run standard NPZD2
import model_NPZD2 as model

days, dt, par, InitCond = model.load_defaults()
output_NPZD2 = model.run(days,dt,InitCond,par)


#%% Run NPZD2_SHELLE coupled model
import model_NPZD2_SHELLE as model

days, dt, par, InitCond = model.load_defaults()
InitCond['Soma'] = 0.
InitCond['Gonad'] = 0.
output_NPZD2_SHELLE = model.run(days,dt,InitCond,par)

#%% Plot
import matplotlib.pyplot as plt

fig, (ax) = plt.subplots(1,1,figsize=(15, 4))
#ax.plot(output_NPZD2['time']/365,output_NPZD2['Phy'],'b-')
#ax.plot(output_NPZD2_SHELLE['time']/365,output_NPZD2_SHELLE['Phy'],'r-')
ax.plot(output_NPZD2['time']/365,output_NPZD2['Zoo'],'b-')
ax.plot(output_NPZD2_SHELLE['time']/365,output_NPZD2_SHELLE['Zoo'],'r-')
ax.set_xlabel('Time (years)')
ax.set_ylabel('Phy (mmol N m$^{-3}$)')
plt.legend(['NPZD2 model','NPZD2_SHELLE coupled model'])