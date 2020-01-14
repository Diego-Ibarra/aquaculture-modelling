# Import Pandas module
import pandas as pd

# Define URL of LOBO buoy
URL = 'http://glider.ceotr.ca/data/live/otn201_sci_water_temp_live.csv'

# Read data from LOBO buoy
glider_data = pd.read_csv(URL,sep=',')

glider_data = glider_data[(glider_data.sci_water_temp < 40) & (glider_data.sci_water_temp > -8)]

glider_data.plot.scatter('unixtime',
                        'depth',
                        c='sci_water_temp',
                        marker='o',
                        edgecolor='none',
                        cmap='jet',
                        ylim=[150,0],
                        figsize=[15,5])