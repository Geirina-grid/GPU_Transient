
# coding: utf-8

# In[1]:

import csv  
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd


# In[2]:

import csv  
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

#The csv files we get from QS is not readable when applying to graphsql,
#So we need to convert it into uft-8 to make it readable.
import sys
reload(sys)
sys.setdefaultencoding('utf-8')

#This function is used for converting csv files which we get from QS into readable csv files. 
def transfer_utf_component(component):
    convert = pd.read_csv(component,encoding="gbk")
    convert.to_csv(component,mode='w',index=False)

    
if __name__ == "__main__":
    transfer_utf_component("Generator.csv")
    transfer_utf_component("Load.csv")
    # transfer_utf_component("Two_winding_transformer.csv")
    # transfer_utf_component("Three_winding_transformer.csv")
#    transfer_utf_component("BusMapping_linux.csv")
#    transfer_utf_component("BaseValue.csv")
#    transfer_utf_component("Breaker.csv")
#    transfer_utf_component("Bus.csv")
#    transfer_utf_component("Compensator_P.csv")
#    transfer_utf_component("Compensator_S.csv")
#     transfer_utf_component("Converter.csv")
#     transfer_utf_component("DCline.csv")
#    transfer_utf_component("Disconnector.csv")
#    transfer_utf_component("Island.csv")
#    transfer_utf_component("Load.csv")
#    transfer_utf_component("Substation.csv")
#     transfer_utf_component("TopoNode.csv")
#    transfer_utf_component("Transformer.csv")
#    transfer_utf_component("Unit.csv")

