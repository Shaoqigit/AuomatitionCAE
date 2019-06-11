# -*- coding: utf-8 -*-
"""
Created on Tue Mar 19 14:42:16 2019

@author: i-shawu
"""

import os
import numpy as np

def lance_abaqus_and_py():

    cmd = "abaqus cae noGUI=Z:\shaoqi-WU\Shaoqi_Stage\modele\Shaoqi_Stage\model_final\Extract_Odb.py"  # the commande will execute in cmd abaqus
    os.system(cmd)
    
#def treat_result()
    
lance_abaqus_and_py()
result = np.loadtxt(open('res.txt','rb'), delimiter=',',skiprows=0)
force_max = np.max(result)

print(force_max)



