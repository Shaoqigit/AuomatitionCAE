# -*- coding: utf-8 -*-
"""
Created on Wed May 02 09:39:53 2018
20/02/2019 modified by shaoqi WU about eveluation points
21/02/2019 modified by shaoqi WU about sigma and calculation of U
25/02/2019 modified by Shaoqi WU about sigma calculated directly in library Pykriging: predict_var
there are different libraries of meta model kriging in internet, we dont know which one? 
27/02/2019 1.stop condition "if" modified by loop while u_min < 2 and add another stop condition : 
            if actual pf equal to the last one, we stop the loop. 
            2. add a comparasion between monte carlo method and AKMCS by calculating error by Shaoqi WU
01/03/2019 training the meta model kriging, should use the real-world values instead units value
07/03/2019 change the way to choose Monte Carlo
"""
# all the libraries we will use are imported here
import numpy as np
from numpy import cosh, sinh

import pyKriging  # library of the Gaussien process
from pyKriging.krige import kriging

import random

import matplotlib.pyplot as plt # library for plotting

from scipy.optimize import minimize, root
import math

import time

#start = time.clock()
# *** Joint properties ***#
# Constants
b = 30.             # Joint width [mm]                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                
l = 60.             # Joint lengthlongueur du joint [mm]
Ep = 62000.         # Young modulus of joined parts [MPa]
ep1 = 4.            # Part 1 thickness [mm]
ep2 = 4             # Part 2 thickness [mm]
F = 6000.          # Applied force [N]
Tau_M = F/(b*l)     # Mean shear stress [MPa]
Tau_R = 27.        # Shear resistance [MPa]

# performance funtion in our problem
def volkersen(s):
    """return the value of performance function which is around 0
    s is the matrix of variable population"""
    X = np.arange(-l/2., l/2.+l/50., l/50.)  # discretize the length of assembly
    
    ea = s[0]
    Ga = s[1]
    
    w = ((Ga/ea)*(1./(Ep*ep1)+1./(Ep*ep2)))**0.5
    Tau = Tau_M *(w*l/2.)*((cosh(w*X)/sinh(w*l/2.))+((Ep*ep1 - Ep*ep2)/(Ep*ep1 + Ep*ep2))*(sinh(w*X)/cosh(w*l/2.)))
    Tau_Max = np.max(Tau)
    return Tau_R - Tau_Max

# have not be used right now, ignore it    
def plot(X, Y):
    minX = int(min(X)/10)*10-10
    maxX = int(max(X)/10)*10+10
    minY = int(min(Y)/10)*10-10
    maxY = int(max(Y)/10)*10+20
    plt.figure(figsize=(10,10))
    plt.plot(X, Y, 'r-o', label='Volkersen',linewidth=2)
    plt.xticks(range(minX,maxX+10,10))
    plt.xlim(minX, maxX)
    plt.xlabel("Distance [mm]",fontsize=22)
    plt.yticks(range(minY,maxY+10,10))
    plt.ylim(minY, maxY)
    plt.ylabel("Shear stress [MPa]",fontsize=20)
    plt.tick_params(axis='both',labelsize=22)
#    plt.legend(loc='upper left',fontsize=18)
#    plt.grid()
        
# sample the variable follow different distribution law
def samplingMC(nMC, style):
    """nMc: number of samples we want to test
       style: it defines the type of the distribution law
       return the samples in determined law"""
    X = []
    if style == 'thickness':
        for x in range(nMC):
            random.seed(a=None)
            x = random.gauss(1, 0.2)
            X.append(x)
    else:
        for x in range(nMC):
            random.seed(a=None)
            x = random.gauss(4000, 800)
            X.append(x)
    
    return X


def plotSamples(S, area, color, marker):
    """print the points"""
    try:
        X1, X2 = [], []
        for s in S:
            X1.append(s[0])
            X2.append(s[1])
        
        scatter = plt.scatter(X1, X2, s = area, c = color, marker=marker)
        plt.xticks(range(-1,2.2,0.1))
        plt.xlim(-1, 2.2)

        plt.yticks(range(3000,5000,100))
        plt.ylim(3000, 5000) 
        plt.xlabel("Thickness",fontsize=15)
        plt.ylabel("Shear modulus",fontsize=15)
    except:
        pass
    return scatter


def plotFunction(X, Y, Line, Label):
    plt.plot(X,Y,Line,label=Label,linewidth=2,markersize=8)    

    
def getInitialPop(S, nMC, N1):
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
            S1.append(S[index])
        if len(Indices) == N1:
            break
       
    S2 = np.delete(S,Indices,axis=0)
    return np.array(S1), np.array(S2), np.array(Indices)
    
# comprend pas pourquoi cette formule de transformation
def convertCtrReal(X, a, b):
    """transfer the normal distribution (0,1) into the real-world value"""
    return a*X+b


def funcFit(X, r1, r2):
    return r1*X+r2

#for interploate the limite state
def funcSolPoly(x):
    maxDeg = len(POLY)
    value = 0
    for i in range(maxDeg):
        deg = float(maxDeg-i-1)
        coef = POLY[i]
        value += coef*x**deg
    return value


def plot_prob(S_POS, S_NEG, S_train=None, S_add=None):
    """plot the probability figures"""

    plt.figure(figsize=(10,10))
    scatter_1 = plotSamples(S_POS, 8, 'lightgreen', 'o')
    scatter_2 = plotSamples(S_NEG, 8, 'gold', 'o')
    plt.legend(handles=[scatter_1, scatter_2],
                   labels=['Positive','Negative'], loc='best')
    if S_train is not None:
        scatter_3 = plotSamples(S_train, 80, 'blue', '*')
        plt.legend(handles=[scatter_1, scatter_2, scatter_3],
                   labels=['Positive', 'Negative', 'training points'], loc='best')
    if S_add is not None:
        scatter_4 = plotSamples(S_add, 150, 'red', '*')
        plt.legend(handles=[scatter_1, scatter_2, scatter_3, scatter_4],
                   labels=['Positive', 'Negative', 'training points', 'Added point'], loc='best')
        
    plt.xlabel("Thickness",fontsize=20)
    plt.ylabel("Shear modulus",fontsize=20)
    plt.tick_params(axis='both',labelsize=15)
        
        

# simulation of Monte Carlo for all points
def RunMcs(S0):
    """simulation of Monte Carlo for all points
    return probability of failure by calculating all the samples"""
    S_NEG = []
    S_POS = []
    Ym = np.zeros(len(S0))
    for enu, s0 in enumerate(S0):     
        Ym[enu] = volkersen(s0)
    
        if Ym[enu] < 0:
            S_NEG.append(s0)
        else:
            S_POS.append(s0)
              
    pf = float(len(S_NEG))/float(len(S0))
    print("failure probability by Monte Carlo = "+str(pf))     
    # plot figures
    plot_prob(S_POS, S_NEG)

    return pf

  
# ------------------------------- Main scipt ---------------------------------# 
# define sampling size
nMC = 3000

# define initial number of individuals to be used to build kriging model 
N1 = 6  # in this case, its enough with 6 points

# sampling nMC value of e0 and G0
E0 = np.array(samplingMC(nMC, 'thickness'))
G0 = np.array(samplingMC(nMC, 'G'))

#plt.scatter(E0, G0)

# convert to real variation domain
#aE, bE = 0.2, 1
#aG, bG = 1000., 4000.  # ?
#E = convertCtrReal(E0,aE,bE)
#G = convertCtrReal(G0,aG,bG)

# reshape inputs, shape(N1,2)
S0 = np.append([E0],[G0],axis=0)
S0 = S0.T
plt.figure(figsize=(10,10))
plotSamples(S0, 10, 'g', 'o')
plt.xlabel("Thickness",fontsize=20)
plt.ylabel("Shear modulus",fontsize=20)
plt.tick_params(axis='both',labelsize=15)

# get initial population (N1 individual)
S1, S2, Indices1 = getInitialPop(S0, nMC, N1)

# run the function of simulation Monte Carlo
S0 = np.array(S0)
t1 = time.time()
pf_Mc = RunMcs(S0)
#print("Pool took : ", time.time() - t1)



# *********************** BUILDING INITIAL KRIGING MODEL *********************#

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
       
    def calcul_prob(S0):
        """calculate the probability of failure by evaluating the sign of evaluation value"""
        S_NEG = []
        S_POS = []
        for s in S0:
                ev_s = k.predict(s)
                if ev_s<=0:
                    S_NEG.append(s)
                else:
                    S_POS.append(s)
        pf = float(len(S_NEG))/float(len(S0))
        print("failure probability= "+str(pf))
        
        return pf, S_POS, S_NEG

    def calcul_u(s):
        """calculate the learning funtion U"""
        ev_s = k.predict(s)
        sigma = math.sqrt(k.predict_var(s))  # k.predict_var si the variance of kriging
        u = abs(ev_s)/sigma
        return sigma, u, ev_s
    
    Y1 = []
    for i in range(N1):
        s1 = S1[i]
        y1 = volkersen(s1)
        Y1.append(y1)
        
    
    k = kriging(S1, Y1, name='simple')
    k.train()
    k.plot()
    ev_matrix = [k.predict(s) for s in S1]
    # ********************************** PLOT ************************************#
    print("\n****  Initialization  ****")
    NEW_S1 = S1
    NEW_S2 = S2
    
    # calculate failure probability 
    pf_matrix = []
    pf, S_POS, S_NEG = calcul_prob(S0)
    pf_matrix.append(pf)
    
        
    plot_prob(S_POS,S_NEG, S1)
    
# ******************************* ITERATION **********************************# 
    u_min = 1.9  #in case of the initial u_min > 5
    iteration = 0
    
    U_min = []
    while u_min < 2.:   
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
            
        u_min = min(U)
        ind_s2min = U.index(u_min)
        G_add = G[ind_s2min]
        ev_matrix.append(G_add)
        U_min.append(u_min)
        print("learning function U : ", u_min)
       
        add_s = S2[ind_s2min]
        NEW_S1 = np.append(S1,[add_s],axis=0)
        print("thickness and moudul: ", add_s)
        print("number of construction ponits: "+str(len(NEW_S1)))
        NEW_S2 = np.delete(S2,ind_s2min,axis=0)
    #   print("add_s= "+str(add_s))
        y1 = volkersen(add_s)
        Y1.append(y1)
        k.addPoint(add_s,y1)
        k.train()
        #k.plot(show=True)
    
    # calculate failure probability
        pf, S_POS, S_NEG = calcul_prob(S0)
        # plot figures
        plot_prob(S_POS, S_NEG, NEW_S1, [add_s])    
        pf_matrix.append(pf)
        
        '''X_plot1 = np.arange(-3., 3.1, 0.1)
        X_plot2 = []
        nPlot = len(X_plot1)
            
        for i in range(nPlot):
            Y_SOL = []
            for j in range(nPlot):
                s_predict = [X_plot1[i],X_plot1[j]]
                y_sol = k.predict(s_predict)
                Y_SOL.append(y_sol)
            POLY = np.polyfit(X_plot1, Y_SOL,3)  # coefficient
            #print(POLY)
            x_sol = root(funcSolPoly, 0.)
    #       print("x_sol= "+str(x_sol.x[0]))
            X_plot2.append(x_sol.x[0])
        plotFunction(X_plot1, X_plot2, 'c-', 'Kriging model')
        if iteration >= 2 and pf_matrix[iteration-1] == pf_matrix[iteration-2]:
            break'''
#    print(ev_matrix)
    return pf_matrix, ev_matrix, U_min
    


pf_matrix, ev_matrix, U_min = AK_MCS(S1, S2, N1, S0)
pf_Mc = RunMcs(S0)

error_matrix = [(pf / pf_Mc) for pf in pf_matrix]
    #print("error between method MC and AKMCS: ",error)
pf = pf_matrix[len(pf_matrix)-1]    
covp = math.sqrt((1 - pf)/(pf * nMC))   
if covp < 0.05:

    print ("C.O.V.P is satified: ", covp)  # variant coefficient of probability
else:
    print("C.O.V.P is : ", covp)


# this part is about plotting the change of the indicator every iteration
def plot_iteration(matrix, style=None):
    """style : g the value of performance function
    e the probability of failure
    return the trend figure plot"""
#print(shape(iteration))
    plt.figure(figsize=(10,10))
    points1 = []
    iterations1 = []
    points2 = []
    iterations2 = []
    
    if style == 'g':        
        for iteration, element in enumerate(matrix):
            if iteration < N1:
                points1 += [element]
                iterations1 += [iteration]
                
            else:
                points2 += [element]
                iterations2 += [iteration]
        point_1, = plt.plot(iterations1, points1, 'b*', markersize=15)
        point_2, = plt.plot(iterations2, points2, 'r*', markersize=15)
        plt.legend(handles=[point_1, point_2,], labels=['Initial training point','Added training point'], loc='best')
        plt.xticks(range(0, len(ev_matrix), 1))
        plt.xlim(-1, len(matrix))
        plt.xlabel('number of training point', fontsize=20)
            #plt.xlim(0, len(ev_matrix))  
        plt.ylabel('value of prediction G(s) at each training point', fontsize=20)
        plt.hlines(0, -1, len(ev_matrix), color = 'c', linestyles = 'dashed')
        
    elif style == 'e':
        for iteration, element in enumerate(matrix):
            if iteration < 1:
                points1 += [element]
                iterations1 += [iteration]
                
            else:
                points2 += [element]
                iterations2 += [iteration]
        point_1, = plt.plot(iterations1, points1, 'b*', markersize=15)
        point_2, = plt.plot(iterations2, points2, 'r*', markersize=15)
        plt.legend(handles=[point_1, point_2],labels=['Initial training point','added training point'], loc='best')
        plt.xticks(range(0, len(matrix), 1))
        plt.xlim(-1, len(matrix))
        plt.xlabel('iteration', fontsize=20)
        plt.ylabel('Normalised Pf',fontsize=20)
        plt.hlines(1, -1, len(matrix), color = 'c', linestyles = 'dashed')
    
    else:
        for iteration, element in enumerate(matrix):
            points1 += [element]
            iterations1 += [iteration]

            point_1, = plt.plot(iterations1, points1, 'b*', markersize=15)
            #point_2, = plt.plot(iterations2, points2, 'r+', markersize=10)
            #plt.legend(handles=[point_1, point_2],labels=['initial training point','added training point'], loc='best')
            plt.xticks(range(0, len(matrix), 1))
            plt.xlim(-1, len(matrix))
            plt.xlabel('iteration', fontsize=20)
            plt.ylabel('value of learning function U ', fontsize=20)
            plt.hlines(2, -1, len(matrix), color = 'c', linestyles = 'dashed')
    #plt.legend(handles=l, labels='limite state', loc='best')
    plt.show()
    
plot_iteration(ev_matrix, 'g')
plot_iteration(error_matrix, 'e')
plot_iteration(U_min)



