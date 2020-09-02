import os
import csv
import numpy as np
import pandas as pd
from collections import defaultdict
import pdb
pd.options.mode.chained_assignment = None  # default='warn'

print('Form Unit Parameter for Case 10790...')

# specify root dir
#   data_epri_SC_wo_hvdc
#root = 'system'
#threshold = 0.000001

# read data

#Unit = pd.read_csv(root + '/Unit.csv')
#UnitPara = pd.read_csv(root + '/Unit_Parameter.csv')
# root_path='g:/Dropbox/summer_intern/IEEE_FILE/'root_path+
print(os.getcwd())
Unit = pd.read_csv('Case10790_Nodeinfo_YT.csv') #raw input Unit.csv
UnitPara = pd.read_csv('Case10790_Unit_Parameter_Scihuan_EPRI.csv') # Unit_Parameter data from EMS sichuan EPRI system renamed to Unit_Parameter_Scihuan_EPRI
#extract only PV-gen bus and slack bus(slack is also generator)
# Unit = Node_info[Node_info['type'] >= 2 ] #or Node_info['type'] == '3'

Unit['Par_Group'] = ''
Unit['Name_Ctrl'] = ''
Unit['Type_Ctrl'] = ''
Unit['Value_Ctrl'] = ''
Unit['Gen_Model'] = ''
Unit['Gen_Par'] = ''
Unit['AVR_Model'] = ''
Unit['AVR_Par'] = ''
Unit['GOV_Model'] = '' 
Unit['GOV_Par'] = '' 
Unit['PSS_Model'] = '' 
Unit['PSS_Par'] = '' 
Unit['Rate_MVA'] = '' 
Unit['Rate_MW'] = '' 
Unit['Xdp'] = '' 
Unit['Xdpp'] = '' 
Unit['X2'] = '' 
Unit['Tj'] = '' 
Unit['P_rate'] = Unit['generation_p']*200.
# Unit.to_csv('Generator_wGen.csv',mode='w', index=False)
# Unit=pd.read_csv('Generator_wGen.csv')
print('Unit Parameter Tabled Formed...')
for iUnit in range(1, len(Unit)):
    #print('Processing Unit no. ' + str(iUnit)+ ' P_rate = ' + str(Unit['P_rate'][iUnit]))
    # pdb.set_trace()
    if Unit['type'][iUnit] >= 2 :
        for iUnitPara in range(1, len(UnitPara)):
            # pdb.set_trace()
            #print('Search Unit Parameter no. ' + str(iUnitPara) + ' Rate_MW = ' + str(UnitPara['Rate_MW'][iUnitPara]))
            if Unit['P_rate'][iUnit] >=UnitPara['Rate_MW'][iUnitPara]: 
                # float(Unit['P_rate'][iUnit]) >=UnitPara['Rate_MW'][iUnitPara]: 
                #print(Unit['Par_Group'][iUnit])
                #print(UnitPara['Par_Group'][iUnitPara])
                Unit['Par_Group'][iUnit] = UnitPara['Par_Group'][iUnitPara]
                Unit['Name_Ctrl'][iUnit] = UnitPara['Name_Ctrl'][iUnitPara]
                Unit['Type_Ctrl'][iUnit] = UnitPara['Type_Ctrl'][iUnitPara]
                Unit['Value_Ctrl'][iUnit] = UnitPara['Value_Ctrl'][iUnitPara]
                Unit['Gen_Model'][iUnit] = UnitPara['Gen_Model'][iUnitPara]
                Unit['Gen_Par'][iUnit] = UnitPara['Gen_Par'][iUnitPara]
                Unit['AVR_Model'][iUnit] = UnitPara['AVR_Model'][iUnitPara]
                Unit['AVR_Par'][iUnit] = UnitPara['AVR_Par'][iUnitPara]
                Unit['GOV_Model'][iUnit] = UnitPara['GOV_Model'][iUnitPara]
                Unit['GOV_Par'][iUnit] = UnitPara['GOV_Par'][iUnitPara]
                Unit['PSS_Model'][iUnit] = UnitPara['PSS_Model'][iUnitPara]
                Unit['PSS_Par'][iUnit] = UnitPara['PSS_Par'][iUnitPara]
                Unit['Rate_MVA'][iUnit] = UnitPara['Rate_MVA'][iUnitPara]
                Unit['Rate_MW'][iUnit] = UnitPara['Rate_MW'][iUnitPara]
                Unit['Xdp'][iUnit] = UnitPara['Xdp'][iUnitPara]
                Unit['Xdpp'][iUnit] = UnitPara['Xdpp'][iUnitPara]
                Unit['X2'][iUnit] = UnitPara['X2'][iUnitPara]
                Unit['Tj'][iUnit] = UnitPara['Tj'][iUnitPara]
                #print('Find the parameter at line ' + str(iUnitPara))
                break;

#new_root = root + '-modified'  
Unit = Unit[Unit['type'] >= 2 ]       
Unit.to_csv('Generator_wGen.csv',mode='w', index=False)
print('Unit Parameter Tabled Formed...')
