# -*- coding: utf-8 -*-
"""
Created on Thu Mar 21 08:53:22 2019

@author: Shaoai WU
This program study the bonded joint with cohesive zone models calculated with
Abaqus .inp
Two variables are modulus and thickness of joint influencing mechanical behavior/
all the calculation resultats is extracted in .odb
"""
import numpy as np
import random
import math

import os, sys

import pyKriging  
from pyKriging.krige import kriging

import scipy
from scipy.optimize import minimize, root


import matplotlib.pyplot as plt

def creat_bat():
    """creat a .bat to run calculation of Abaqus without opening CAE
    we can use it to run paramellely many jobs"""
    
    # the command is written in .bat
    command = "@echo off\n\
    cd G:\CODE\code_stage\modele_final\n\
    abq2016 job=NewInp cpus=4 scratch=G:\CODE\code_stage\model_final\ interactive ask_delete=OFF"
    bat_name = main_dir + "batFile.bat"
    bat_file = open(bat_name, 'w')
    bat_file.write(command)
    bat_file.close()
    
    #delete those space in front of every line in the file.bat
    bat_file = open(bat_name, 'r+')
    lines = bat_file.readlines()
    bat_file.close()
    bat_file = open(bat_name, 'w')
    bat_file.writelines(line.lstrip() for line in lines)
    bat_file.close()

# functions for plotting figure
def plotSamples(S, area, color):
    """print the points"""
    try:
        X1, X2 = [], []
        for s in S:
            X1.append(s[0])
            X2.append(s[1])
        
        plt.scatter(X1, X2, s = area, c = color)
        plt.xticks(range(-1,2.2,0.1))
        plt.xlim(-1, 2.2)

        plt.yticks(range(3000,5000,100))
        plt.ylim(3000, 5000) 
    except:
        pass

def plot_prob(S_POS, S_NEG, S_train=None, S_add=None):
    """plot the probability figures"""

    plt.figure(figsize=(10,10))
    plotSamples(S_POS, 10, 'g')
    plotSamples(S_NEG, 10, 'r')
    if S_train is not None:
        plotSamples(S_train, 40, 'bs')
    if S_add is not None:
        plotSamples(S_add, 80, 'y')
    plt.xlabel("Thickness",fontsize=20)
    plt.ylabel("Shear modulus",fontsize=20)
    plt.tick_params(axis='both',labelsize=15)


def sampling_MC(nMC, style):
    """nMc: number of samples we want to test
       style: it defines the type of the distribution law
       return the samples in determined law"""
    X = []
    if style == 'thickness':
        for x in range(nMC):
            random.seed(a=None)
            x = random.gauss(0.5, 0.1)
            X.append(x)
    else:
        for x in range(nMC):
            random.seed(a=None)
            x = random.gauss(Ec, 800)
            X.append(x)
    
    return X
 
    
def get_initial_pop(S0, nMC, N1):
    """selection the first population in Monte Carlo
       return S1: N1 points in it
       S2 : nMC - N1 point in it
       indices: index"""
    Indices = []
    S1, S2 = [], []
    for i in range(1000*nMC):
        random.seed(a=None)
        index = random.randint(0,nMC)
        if index not in Indices and len(Indices)<N1:
            Indices.append(index)
            S1.append(S0[index])
        if len(Indices) == N1:
            break
       
    S2 = np.delete(S0,Indices,axis=0)
    return np.array(S1), np.array(S2), np.array(Indices)

def coord_node(node):
    if node - int(node) == 0:
        coord = str(int(node))+'.'
    else:
        coord = str(node)
    return coord

def format_mc(S):
    """formate the population of thickness and Modulu Young in file .csv"""
    texte = ''
    for ii_A in range(0,len(S[:,0])):
       t,E = coord_node(S[ii_A,0]),coord_node(S[ii_A,1])
       texte +='{0:>3s}'.format(t)+','+ '{0:>10s}'.format(E) +  '\n'
    return texte

def cohesive_stiffness(t, E):
    """formule in these LE_GOFF page 146"""
    Knn = E/t  # Knn stiffness of cohesive interface in traction
    G = E/(2*(1+v))  # transversal modulu
    Kss, Ktt = G/t, G/t  # shear stiffness in cohesive interface
    return int(Knn), int(Kss), int(Ktt)


def damage_init(t):
    stress_init_1 = -14.4 * t + 41.58 # approximate from data in article
    stress_init_2 = stress_init_1 * 0.95
    return int(stress_init_1), int(stress_init_2)

def damage_energy(t):
    energy_23 = 780.8 * t + 406.6 # idem damage_energy
    energy_1 = energy_23 / 2.5
    return energy_1, energy_23
    
def replace_text(old_text, new_text, old_inp_text):
    """replace the content in text of file.inp"""
    new_inp_text = old_inp_text.replace(old_text, new_text)
    return new_inp_text


def modify_inp(ts, Es):
    """modify the stiffness with creating thickness and modulu Young in file.inp"""
    old_inp_name = main_dir + "czm_a_m40.inp"  #search old .inp file
    old_inp = open(old_inp_name,'r')
    old_inp_text = old_inp.read()  # f.read([size]) return a list of [size] dimension in file
    
    new_inp_name = main_dir + "NewInp.inp"
    new_inp = open(new_inp_name, 'w')
    
    old_text = "8000.,2800.,2800."
    Knn, Kss, Ktt = cohesive_stiffness(ts, Es)
    new_text = str(Knn) + "," + str(Kss) + "," + str(Ktt)
    print("stiffness of cohesive interface:\nKnn   Kss   Ktt\n", Knn, Kss, Ktt)
    
    new_inp_text = replace_text(old_text, new_text, old_inp_text)
    new_inp.write(new_inp_text)
    
    old_inp.close()
    new_inp.close()
    

def lanch_calcul(ts, Es):
    """This function is very important like the volkersen function in script AKMCS_f2_MC _final.py
    as a performance function.
    here criterion is maximum force of failure, return difference between resistance of f_m and compute f_m
    launch a calculation of Abaqus without opening windows CAE by using file.bat
    and extract the calculation results in file.odb"""
    #modify_inp(ts, Es)  # obtain the modified .inp with new stiffness
    
    
    os.chdir(main_dir)
    os.system(main_dir + "batFile.bat")   
    
    # lanch the program python for extract data in .odb without opening Abaqus CAE
    cmd = "abaqus cae noGUI=Z:\shaoqi-WU\Shaoqi_Stage\modele\Shaoqi_Stage\model_final\Extract_Odb.py"
    os.system(cmd)
    
    result = np.loadtxt(open("res.txt",'rb'), delimiter=',', skiprows=0)
    force_max = np.max(result)
    
    Y = float(f_r - force_max)
    return Y
    
def write_S(S, name):
    """functon for writing the matrix of population in file.txt"""
    pop_S = main_dir + name + ".csv"
    S_csv = open(pop_S, 'w')
    pop_S_text = format_mc(S)
    S_csv.write(pop_S_text)
    S_csv.close()
    

#----------------------------------------------main script-----------------------------------#
# constant of meta-model
nMC = 3000  # number of Monte Carlo
N1 = 8  # number of initial traning points

# path of all the files saved
main_dir = "Z:\\shaoqi-WU\\Shaoqi_Stage\\modele\\Shaoqi_Stage\\model_final\\"  # set actual working path
os.chdir(main_dir)  # enter actual working path
#print(os.getcwd())  #check the actual working parth

# parameters of ahesive
Ec = 4000  # Young Modulu in MPa
v = 0.49  # poisson coefficient 
Rc = 21  # elastic limitation
f_r = 8.90  # force ressistance

# initialization the sampling of varialbes
t = np.array(sampling_MC(nMC,'thickness'))  # varialebs of adhesive thickness
E = np.array(sampling_MC(nMC, 'module young'))  # variable of adhsive modulu Yong
 
# reshape in dimension (nMC, 2)
S0 = np.append([t],[E],axis=0)
S0 = S0.T  # transposition

# get initial population (N1 inividual)
S1, S2, Index = get_initial_pop(S0, nMC, N1)

# create files to saving populations
write_S(S0, "S0")
write_S(S1, "S1")
write_S(S2, "S2")

#g = lanch_calcul(main_dir, 0.7, 4200)

def AK_MCS(S1, S2, N1, S0):
    """main function for AK-MCS: active learning method based on meta-model kriging
    with simulation of Monte Carlo
    S1: first selection population in Monte Carlo to construct the meta-model
    S2: explained before
    S0: the whole population generated by Monte Carlo
    return 
    pf_matrix: the list of all calculated probabilities in every loop
    ev_matrix: the list of all evaluation value of performance function by meta-model
    U_min: the list of all minimum learning functions U in every loop
    return them for the data treatment"""
    # INITIALIZATION
    #constraction of first model of Kriging
    
    def calcul_prob_kriging(S0):        
        """calculate the probability of failure by evaluating the sign of evaluation value"""
        S_NEG = []
        S_POS = []
        for s in S0:
                ev_s = k.predict(s)
                if ev_s<=0:
                    S_NEG.append(s)
                else:
                    S_POS.append(s)
        pf = round(float(len(S_NEG))/float(len(S0)),4)  # take 4 number behind the .
        print("failure probability= "+str(pf))
        
        return pf, S_POS, S_NEG
    def calcul_u(s):
        """calculate the learning function U"""
        pred_s = k.predict(s)
        sigma = math.sqrt(k.predict_var(s))  # k.predict_var is the variance of Kriging
        u = abs(pred_s)/sigma
        return sigma, u, pred_s
    
    
     # results of first calcultion
    t1 = S1[:, 0]  # thickness in S1
    E1 = S1[:, 1]  # modulu Young in S1
    fichier = open('result.txt','w')
    fichier.close()  
    
    for ts, Es in zip(t1, E1):
        Y = lanch_calcul(ts, Es)  
        
        result = open("result.txt",'a')  # 'a' means open a file for append(write only str) 
        result.write(str(Y)+'\n')
        result.close()
        
    result_name=main_dir + "result.txt"
    result_file=open(result_name,"r")
    Y=np.loadtxt('result.txt')  # 
    result_file.close()    
    # library of kriging
    k = kriging(S1, Y, name='simple')
    k.train()
    pred_matrix = [k.predict(s) for s in S1]
    
    #first calculation of probability
    print("\n******* Initialization********")
    # calculate failure probability
    pf, S_POS, S_NEG = calcul_prob_kriging(S0)
    plot_prob(S_POS, S_NEG, S1)

#AK_MCS(S1, S2, N1, S0)
    #ITERATION
    u_min = 1.9  #in case of the initial u_min > 5
    iteration = 0
    pf_matrix = []
    U_min = []
    NEW_S1 = S1
    NEW_S2 = S2
    while u_min < 2:
        S1 = NEW_S1
        S2 = NEW_S2
        U = []
        G = []
        iteration += 1
        print("\n****  iteration number: "+str(iteration)+ "****")
        
        for j in range(0, len(S2), 1):
            s2 = S2[j]    
            sigma, u, ev_s = calcul_u(s2)
            U.append(u)
            G.append(ev_s)
            
        u_min = min(U)  # search for the min U
        ind_s2min = U.index(u_min)
        G_add = G[ind_s2min]
        pred_matrix.append(G_add)
        U_min.append(u_min)  # collect all the Umin in every interation
        print("learning function U : ", u_min)
        
        add_s = S2[ind_s2min]
        NEW_S1 = np.append(S1, [add_s], axis=0)
        NEW_S2 = np.delete(S2, ind_s2min, axis=0)
        
        
        t_add = add_s[0]
        E_add = add_s[1]
        add_Y = lanch_calcul(t_add, E_add)
        k.addPoint(add_s, add_Y)
        k.train()
        
    #  calculation failure probability
        pf, S_POS, S_NEG = calcul_prob_kriging(S0)
        if (iteration % 10) == 0:
            plot_prob(S_POS, S_NEG, S1)
        else:
            pass
        
        pf_matrix.append(pf)
    return pf_matrix, U_min
        
            
pf_matrix, U_min = AK_MCS(S1, S2, N1, S0)        
    




































