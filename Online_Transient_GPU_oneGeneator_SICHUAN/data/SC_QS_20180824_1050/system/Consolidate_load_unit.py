import os
import csv
import numpy as np
import pandas as pd
from collections import defaultdict
pd.options.mode.chained_assignment = None  # default='warn'

print('Consolidata Unit and Load...')

# specify root dir
#   data_epri_SC_wo_hvdc
#   data_epri_36-test
# root = 'result_sc'

#read data
# generator = pd.read_csv(root + '/Unit.csv')
# load = pd.read_csv(root + '/Load.csv')

generator = pd.read_csv('Unit.csv')
load = pd.read_csv('Load.csv')

# Consolidata loads at the same bus
for i in range(1, len(load)-1):
#for i in range(1, 100):
    # print('Processing load =' + str(i))
    if (load['off'][i] == "0" and load['island'][i] != "-1"):
        for j in range(i+1, len(load)):
            if load['bs'][j]==load['bs'][i] and load['off'][j] == "0":
                print('load i=' + str(i) + ' P= ' + str(load['P'][i]) + ' load j=' + str(j) + ' P= ' + str(load['P'][j]))
                load['P'][i] =  str(float(load['P'][i]) + float(load['P'][j])) 
                print('load i=' + str(i) + ' P= ' + str(load['P'][i]) + ' load j=' + str(j) + ' P= ' + str(load['P'][j]))
                load['Q'][i] =  str(float(load['Q'][i]) + float(load['Q'][j])) 
                if (load['P_meas'][i] != "''" and load['P_meas'][j] !=  "''"):
                    load['P_meas'][i] =  str(float(load['P_meas'][i]) + float(load['P_meas'][j])) 
                if (load['Q_meas'][i] !=  "''" and load['Q_meas'][j] != "''"):    
                    load['Q_meas'][i] =  str(float(load['Q_meas'][i]) + float(load['Q_meas'][j]))
                load['off'][j] = "1"
 
# Consolidata unit at the same bus
for i in range(1, len(generator)-1):
#for i in range(40, 50):
    # print('Processing generator =' + str(i))
    #print('1 Generator i=' + str(i) + ' off ' + str(generator['off'][i]) + ' island = ' + str(generator['island'][i]))
    if (generator['off'][i] == "0" and generator['island'][i] != "-1"):
        #print('2 Generator i=' + str(i) + ' off ' + str(generator['off'][i]) + ' island = ' + str(generator['island'][i]))
        for j in range(i+1, len(generator)):
            #print('Generator i=' + str(i) + 'at ' + str(generator['bs'][i]) + 'generator j=' + str(j) + ' at ' + str(generator['bs'][j]))
            if generator['bs'][j]==generator['bs'][i] and generator['off'][j] == "0":
                generator['P'][i] = str(float(generator['P'][i]) + float(generator['P'][j]))
                generator['Q'][i] = str(float(generator['Q'][i]) + float(generator['Q'][j])) 
                generator['P_max'][i] = str(float(generator['P_max'][i]) + float(generator['P_max'][j])) 
                generator['P_min'][i] = str(float(generator['P_min'][i]) + float(generator['P_min'][j])) 
                generator['Q_max'][i] = str(float(generator['Q_max'][i]) + float(generator['Q_max'][j])) 
                generator['Q_min'][i] = str(float(generator['Q_min'][i]) + float(generator['Q_min'][j])) 
                print('Generator i=' + str(i) + ' P_rate= ' + str(generator['P_rate'][i]) + ' generator j=' + str(j) + ' P_rate= ' + str(generator['P_rate'][j]))
                generator['P_rate'][i] = str(float(generator['P_rate'][i]) + float(generator['P_rate'][j]))
                print('After Generator i=' + str(i) + ' P_rate= ' + (generator['P_rate'][i]))
                if (generator['P_meas'][i] != "''" and generator['P_meas'][j] !=  "''"):
                    generator['P_meas'][i] = str(float(generator['P_meas'][i]) + float(generator['P_meas'][j]))
                if (generator['Q_meas'][i] !=  "''" and generator['Q_meas'][j] != "''"):
                    generator['Q_meas'][i] = str(float(generator['Q_meas'][i]) + float(generator['Q_meas'][j]))
                generator['off'][j] = "1"
 
# output data
# new_root = root + '-modified/'
# if not os.path.exists(new_root):
#     os.makedirs(new_root)

# generator.to_csv(new_root + '/Unit.csv',mode='w', index=False)
# load.to_csv(new_root + '/Load.csv',mode='w', index=False)

generator.to_csv('Unit_Consolidated.csv',mode='w', index=False)
load.to_csv('Load_Consolidated.csv',mode='w', index=False)
print(' Unit and Load are consolidated...')

