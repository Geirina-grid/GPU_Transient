# coding: utf-8

# In[2]:

import csv  
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd


# # Two_port transformer

# In[2]:

df = pd.read_csv('Transformer.csv')
df = df.drop([0])


# In[3]:

two_node_transformer = df[df['type'] == '2']


# In[6]:

two_node_transformer['Jtap_V'] = pd.to_numeric(two_node_transformer['Jtap_V'])


# In[7]:

two_node_transformer['Jtap_V'] = two_node_transformer['Jtap_V'].replace(10.5, 0)
two_node_transformer['Jtap_V'] = two_node_transformer['Jtap_V'].replace(6.3, 0)
two_node_transformer['Jtap_V'] = two_node_transformer['Jtap_V'].replace(6, 0)
two_node_transformer['Jtap_V'] = two_node_transformer['Jtap_V'].replace(35, 0)
two_node_transformer['Jtap_V'] = two_node_transformer['Jtap_V'].replace(31.5, 0)
two_node_transformer['Jtap_V'] = two_node_transformer['Jtap_V'].replace(10, 0)
two_node_transformer['Jtap_V'] = two_node_transformer['Jtap_V'].replace(37, 0)
two_node_transformer['Jtap_V'] = two_node_transformer['Jtap_V'].replace(20, 0)


# In[8]:

temp = two_node_transformer.loc[:,'id']
temp = temp.apply(lambda x: pd.Series([float(x+'.1'),float(x+'.2')]))
temp = pd.DataFrame(temp.values.reshape(-1,1,order='C'))
temp.columns = ['index']

two_node_transformer.to_csv('two_port_transformer_Cinput.csv',mode='w',index=False)

three_node_transformer = df[df['type'] == '3']
three_node_transformer['Jtap_V'] = pd.to_numeric(three_node_transformer['Jtap_V'])

three_node_transformer['Jtap_V'] = three_node_transformer['Jtap_V'].replace(10.5, 0)
three_node_transformer['Jtap_V'] = three_node_transformer['Jtap_V'].replace(11, 0)
three_node_transformer['Jtap_V'] = three_node_transformer['Jtap_V'].replace(6, 0)
three_node_transformer['Jtap_V'] = three_node_transformer['Jtap_V'].replace(35, 0)
three_node_transformer.to_csv('three_port_transformer_Cinput.csv',mode='w',index=False)



# In[9]:

name_for_two_transformer = two_node_transformer.iloc[np.repeat(np.arange(len(two_node_transformer)), 2)]
name_two_transformer = name_for_two_transformer[['name']]
name_two_transformer = name_two_transformer.reset_index(drop=False)


# In[10]:

type_for_two_transformer = two_node_transformer.iloc[np.repeat(np.arange(len(two_node_transformer)), 2)]
type_two_transformer = type_for_two_transformer[['type']]
type_two_transformer = type_two_transformer.reset_index(drop=False)


# In[11]:

port_point_volt = two_node_transformer.loc[:,['I_Vol','J_Vol']].copy().values
port_point_S = two_node_transformer.loc[:,['I_S','J_S']].copy().values
port_point_itapH = two_node_transformer.loc[:,['Itap_H','Jtap_V']].copy().values
port_point_itapL = two_node_transformer.loc[:,['Itap_L','Jtap_V']].copy().values
port_point_itapC = two_node_transformer.loc[:,['Itap_C','Jtap_V']].copy().values
port_point_R = two_node_transformer.loc[:,['Ri','Rj']].copy().values
port_point_X = two_node_transformer.loc[:,['Xi','Xj']].copy().values
port_point_node = two_node_transformer.loc[:,['I_node','J_node']].copy().values
port_point_P = two_node_transformer.loc[:,['I_P','J_P']].copy().values
port_point_Q = two_node_transformer.loc[:,['I_Q','J_Q']].copy().values
port_point_off = two_node_transformer.loc[:,['I_off','J_off']].copy().values
# port_point_pimeas = two_node_transformer.loc[:,['Pi_meas','Pj_meas']].copy().values
# port_point_qimeas = two_node_transformer.loc[:,['Qi_meas','Qj_meas']].copy().values
port_point_Tx = two_node_transformer.loc[:,['I_nd','J_nd']].copy().values
port_point_bs = two_node_transformer.loc[:,['I_bs','J_bs']].copy().values
port_point_island = two_node_transformer.loc[:,['I_island','J_island']].copy().values
port_point_Rstar = two_node_transformer.loc[:,['Ri*','Rj*']].copy().values
port_point_Xstar = two_node_transformer.loc[:,['Xi*','Xj*']].copy().values
port_point_t = two_node_transformer.loc[:,['I_t','J_t']].copy().values
port_point_base = two_node_transformer.loc[:,['IBase_V','JBase_V']].copy().values


# In[12]:

try_add_volt = pd.DataFrame(port_point_volt.reshape(-1,1,order='C'))
try_add_S = pd.DataFrame(port_point_S.reshape(-1,1,order='C'))
try_add_itapH = pd.DataFrame(port_point_itapH.reshape(-1,1,order='C'))
try_add_itapL = pd.DataFrame(port_point_itapL.reshape(-1,1,order='C'))
try_add_itapC = pd.DataFrame(port_point_itapC.reshape(-1,1,order='C'))
try_add_R = pd.DataFrame(port_point_R.reshape(-1,1,order='C'))
try_add_X = pd.DataFrame(port_point_X.reshape(-1,1,order='C'))
try_add_node = pd.DataFrame(port_point_node.reshape(-1,1,order='C'))
try_add_P = pd.DataFrame(port_point_P.reshape(-1,1,order='C'))
try_add_Q = pd.DataFrame(port_point_Q.reshape(-1,1,order='C'))
try_add_off = pd.DataFrame(port_point_off.reshape(-1,1,order='C'))
# try_add_pimeas = pd.DataFrame(port_point_pimeas.reshape(-1,1,order='C'))
# try_add_qimeas = pd.DataFrame(port_point_qimeas.reshape(-1,1,order='C'))
try_add_Tx = pd.DataFrame(port_point_Tx.reshape(-1,1,order='C'))
try_add_bs = pd.DataFrame(port_point_bs.reshape(-1,1,order='C'))
try_add_island = pd.DataFrame(port_point_island.reshape(-1,1,order='C'))
try_add_Rstar = pd.DataFrame(port_point_Rstar.reshape(-1,1,order='C'))
try_add_Xstar = pd.DataFrame(port_point_Xstar.reshape(-1,1,order='C'))
try_add_t = pd.DataFrame(port_point_t.reshape(-1,1,order='C'))
try_add_base = pd.DataFrame(port_point_base.reshape(-1,1,order='C'))


# In[13]:

try_add_volt.columns = ['volt']
try_add_S.columns = ['S']
try_add_itapH.columns = ['itapH']
try_add_itapL.columns = ['itapL']
try_add_itapC.columns = ['itapC']
try_add_R.columns = ['R']
try_add_X.columns = ['X']
try_add_node.columns = ['node']
try_add_P.columns = ['P']
try_add_Q.columns = ['Q']
try_add_off.columns = ['off']
# try_add_pimeas.columns = ['pimeas']
# try_add_qimeas.columns = ['qimeas']
try_add_Tx.columns = ['Tx']
try_add_bs.columns = ['bs']
try_add_island.columns = ['island']
try_add_Rstar.columns = ['Rstar']
try_add_Xstar.columns = ['Xstar']
try_add_t.columns = ['t']
try_add_base.columns = ['base']


# In[14]:

temp[['name']] = name_two_transformer[['name']]
temp[['type']] = type_two_transformer[['type']]
temp[['volt']] = try_add_volt[['volt']]
temp[['S']] = try_add_S[['S']]
temp[['itapH']] = try_add_itapH[['itapH']]
temp[['itapL']] = try_add_itapL[['itapL']]
temp[['itapC']] = try_add_itapC[['itapC']]
temp[['R']] = try_add_R[['R']]
temp[['X']] = try_add_X[['X']]
temp[['node']] = try_add_node[['node']]
temp[['P']] = try_add_P[['P']]
temp[['Q']] = try_add_Q[['Q']]
temp[['off']] = try_add_off[['off']]
# temp[['pimeas']] = try_add_pimeas[['pimeas']]
# temp[['qimeas']] = try_add_qimeas[['qimeas']]
temp[['nd']] = try_add_Tx[['Tx']]
temp[['bs']] = try_add_bs[['bs']]
temp[['island']] = try_add_island[['island']]
temp[['Rstar']] = try_add_Rstar[['Rstar']]
temp[['Xstar']] = try_add_Xstar[['Xstar']]
temp[['t']] = try_add_t[['t']]
temp[['base']] = try_add_base[['base']]


# In[15]:

temp.to_csv('two_port_transformer_Cinput.csv',mode='w',index=False)


# # three_port transformer

# In[3]:

df = pd.read_csv('Transformer.csv')
df = df.drop([0])
three_node_transformer = df[df['type'] == '3']


# In[4]:

# three_node_transformer = three_node_transformer[['id','name','type','I_Vol','K_Vol','J_Vol','I_S','K_S','J_S','Itap_H','Itap_L','Itap_C','Ktap_H','Ktap_L','Ktap_C','Jtap_V','Ri','Xi','Rk','Xk','Rj','Xj','I_node','K_node','J_node','I_P','I_Q','K_P','K_Q','J_P','J_Q','I_off','K_off','J_off','Pi_meas','Qi_meas','Pk_meas','Qk_meas','Pj_meas','Qj_meas','I_nd','K_nd','J_nd','I_bs','K_bs','J_bs','I_island','K_island','J_island','Ri*','Xi*','Rk*','Xk*','Rj*','Xj*','I_t','K_t','J_t','IBase_V','KBase_V','JBase_V']]


# In[5]:

temp1 = three_node_transformer.loc[:,'id']
temp1 = temp1.apply(lambda x: pd.Series([float(x+'.1'),float(x+'.2'),float(x+'.3')]))
temp1 = pd.DataFrame(temp1.values.reshape(-1,1,order='C'))

name_for_three_transformer = three_node_transformer.iloc[np.repeat(np.arange(len(three_node_transformer)), 3)]
name_three_transformer = name_for_three_transformer[['name']]
name_three_transformer = name_three_transformer.reset_index(drop=False)

type_for_three_transformer = three_node_transformer.iloc[np.repeat(np.arange(len(three_node_transformer)), 3)]
type_three_transformer = type_for_three_transformer[['type']]
type_three_transformer = type_three_transformer.reset_index(drop=False)


# In[6]:

three_node_transformer['Jtap_V'] = pd.to_numeric(three_node_transformer['Jtap_V'])

three_node_transformer['Jtap_V'] = three_node_transformer['Jtap_V'].replace(10.5, 0)
three_node_transformer['Jtap_V'] = three_node_transformer['Jtap_V'].replace(11, 0)
three_node_transformer['Jtap_V'] = three_node_transformer['Jtap_V'].replace(6, 0)
three_node_transformer['Jtap_V'] = three_node_transformer['Jtap_V'].replace(35, 0)


# In[7]:

port_point_three_volt = three_node_transformer.loc[:,['I_Vol','K_Vol','J_Vol']].copy().values
port_point_three_S = three_node_transformer.loc[:,['I_S','K_S','J_S']].copy().values
port_point_three_itapH = three_node_transformer.loc[:,['Itap_H','Ktap_H','Jtap_V']].copy().values
port_point_three_itapL = three_node_transformer.loc[:,['Itap_L','Ktap_L','Jtap_V']].copy().values
port_point_three_itapC = three_node_transformer.loc[:,['Itap_C','Ktap_C','Jtap_V']].copy().values
port_point_three_R = three_node_transformer.loc[:,['Ri','Rk','Rj']].copy().values
port_point_three_X = three_node_transformer.loc[:,['Xi','Xk','Xj']].copy().values
port_point_three_node = three_node_transformer.loc[:,['I_node','K_node','J_node']].copy().values
port_point_three_P = three_node_transformer.loc[:,['I_P','K_P','J_P']].copy().values
port_point_three_Q = three_node_transformer.loc[:,['I_Q','K_Q','J_Q']].copy().values
port_point_three_off = three_node_transformer.loc[:,['I_off','K_off','J_off']].copy().values
# port_point_three_pimeas = three_node_transformer.loc[:,['Pi_meas','Pk_meas','Pj_meas']].copy().values
# port_point_three_qimeas = three_node_transformer.loc[:,['Qi_meas','Qk_meas','Qj_meas']].copy().values
port_point_three_Tx = three_node_transformer.loc[:,['I_nd','K_nd','J_nd']].copy().values
port_point_three_bs = three_node_transformer.loc[:,['I_bs','K_bs','J_bs']].copy().values
port_point_three_island = three_node_transformer.loc[:,['I_island','K_island','J_island']].copy().values
port_point_three_Rstar = three_node_transformer.loc[:,['Ri*','Rk*','Rj*']].copy().values
port_point_three_Xstar = three_node_transformer.loc[:,['Xi*','Xk*','Xj*']].copy().values
port_point_three_t = three_node_transformer.loc[:,['I_t','K_t','J_t']].copy().values
port_point_three_base = three_node_transformer.loc[:,['IBase_V','KBase_V','JBase_V']].copy().values


# In[8]:

try_add_three_volt = pd.DataFrame(port_point_three_volt.reshape(-1,1,order='C'))
try_add_three_S = pd.DataFrame(port_point_three_S.reshape(-1,1,order='C'))
try_add_three_itapH = pd.DataFrame(port_point_three_itapH.reshape(-1,1,order='C'))
try_add_three_itapL = pd.DataFrame(port_point_three_itapL.reshape(-1,1,order='C'))
try_add_three_itapC = pd.DataFrame(port_point_three_itapC.reshape(-1,1,order='C'))
try_add_three_R = pd.DataFrame(port_point_three_R.reshape(-1,1,order='C'))
try_add_three_X = pd.DataFrame(port_point_three_X.reshape(-1,1,order='C'))
try_add_three_node = pd.DataFrame(port_point_three_node.reshape(-1,1,order='C'))
try_add_three_P = pd.DataFrame(port_point_three_P.reshape(-1,1,order='C'))
try_add_three_Q = pd.DataFrame(port_point_three_Q.reshape(-1,1,order='C'))
try_add_three_off = pd.DataFrame(port_point_three_off.reshape(-1,1,order='C'))
# try_add_three_pimeas = pd.DataFrame(port_point_three_pimeas.reshape(-1,1,order='C'))
# try_add_three_qimeas = pd.DataFrame(port_point_three_qimeas.reshape(-1,1,order='C'))
try_add_three_Tx = pd.DataFrame(port_point_three_Tx.reshape(-1,1,order='C'))
try_add_three_bs = pd.DataFrame(port_point_three_bs.reshape(-1,1,order='C'))
try_add_three_island = pd.DataFrame(port_point_three_island.reshape(-1,1,order='C'))
try_add_three_Rstar = pd.DataFrame(port_point_three_Rstar.reshape(-1,1,order='C'))
try_add_three_Xstar = pd.DataFrame(port_point_three_Xstar.reshape(-1,1,order='C'))
try_add_three_t = pd.DataFrame(port_point_three_t.reshape(-1,1,order='C'))
try_add_three_base = pd.DataFrame(port_point_three_base.reshape(-1,1,order='C'))


# In[9]:

try_add_three_volt.columns = ['volt']
try_add_three_S.columns = ['S']
try_add_three_itapH.columns = ['itapH']
try_add_three_itapL.columns = ['itapL']
try_add_three_itapC.columns = ['itapC']
try_add_three_R.columns = ['R']
try_add_three_X.columns = ['X']
try_add_three_node.columns = ['node']
try_add_three_P.columns = ['P']
try_add_three_Q.columns = ['Q']
try_add_three_off.columns = ['off']
# try_add_three_pimeas.columns = ['pimeas']
# try_add_three_qimeas.columns = ['qimeas']
try_add_three_Tx.columns = ['Tx']
try_add_three_bs.columns = ['bs']
try_add_three_island.columns = ['island']
try_add_three_Rstar.columns = ['Rstar']
try_add_three_Xstar.columns = ['Xstar']
try_add_three_t.columns = ['t']
try_add_three_base.columns = ['base']


# In[10]:

temp1.columns = ['index']
temp1[['name']] = name_three_transformer[['name']]
temp1[['type']] = type_three_transformer[['type']]
temp1[['volt']] = try_add_three_volt[['volt']]
temp1[['S']] = try_add_three_S[['S']]
temp1[['itapH']] = try_add_three_itapH[['itapH']]
temp1[['itapL']] = try_add_three_itapL[['itapL']]
temp1[['itapC']] = try_add_three_itapC[['itapC']]
temp1[['R']] = try_add_three_R[['R']]
temp1[['X']] = try_add_three_X[['X']]
temp1[['node']] = try_add_three_node[['node']]
temp1[['P']] = try_add_three_P[['P']]
temp1[['Q']] = try_add_three_Q[['Q']]
temp1[['off']] = try_add_three_off[['off']]
# temp1[['pimeas']] = try_add_three_pimeas[['pimeas']]
# temp1[['qimeas']] = try_add_three_qimeas[['qimeas']]
temp1[['nd']] = try_add_three_Tx[['Tx']]
temp1[['bs']] = try_add_three_bs[['bs']]
temp1[['island']] = try_add_three_island[['island']]
temp1[['Rstar']] = try_add_three_Rstar[['Rstar']]
temp1[['Xstar']] = try_add_three_Xstar[['Xstar']]
temp1[['t']] = try_add_three_t[['t']]
temp1[['base']] = try_add_three_base[['base']]
temp1.to_csv('three_port_transformer_Cinput.csv',mode='w',index=False)


# # neutral point

# In[12]:

# build neutral point
zhongxingdian = three_node_transformer[['id']]
zhongxingdian = zhongxingdian.reset_index(drop=False)
zhongxingdian['index']*=-1
zhongxingdian.columns = ['middle_point','id']
final_zhongxingdian = zhongxingdian[['middle_point']]
final_zhongxingdian.to_csv('Neutral_point_transformer_Cinput.csv',mode='w',index=False)


# # change_transformer_to_nonzero

# In[13]:

transformer_two = pd.read_csv('two_port_transformer_Cinput.csv')
transformer_two['Xstar'] = pd.to_numeric(transformer_two['Xstar'])
transformer_two['Xstar'] = transformer_two['Xstar'].fillna(0.001)
transformer_two['Xstar'] = transformer_two['Xstar'].replace(0, 0.001)
transformer_two.to_csv('two_port_transformer_Cinput.csv',mode='w',index=False)
transformer_three = pd.read_csv('three_port_transformer_Cinput.csv')
transformer_three['Xstar'] = transformer_three['Xstar'].fillna(0.001)
transformer_three['Xstar'] = transformer_three['Xstar'].replace(0, 0.001)
transformer_three.to_csv('three_port_transformer_Cinput.csv',mode='w',index=False)

