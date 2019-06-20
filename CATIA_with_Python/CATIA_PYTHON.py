# -*- coding: utf-8 -*-
"""
Created on Fri Apr 19 09:58:21 2019

@author: Shaoqi WU
"""
import os
from win32com.client import Dispatch

import numpy as np
import matplotlib.pyplot as plt


def NACA_Profil():
    NACA_type = "2412"
    
    M_init = float(NACA_type[0])
    P_init = float(NACA_type[1])
    T_init = float(NACA_type[2:4])
    
    dis_points = 500
    
    c = 100
    a0 = 0.2969
    a1 = -0.126
    a2 = -0.3516
    a3 = 0.2843
    a4 = -0.1015
    
    M = M_init/100
    P = P_init/10
    T = T_init/100
    
    
    yc = np.ones(shape=(500, 1))
    dyc_dx = np.ones(shape=(500, 1))
    theta = np.ones(shape=(500, 1))
    
    x = np.linspace(0, 1, dis_points)
    for i in range(500):
        if x[i] >= 0 and x[i] < P:
            yc[i] = (M/P**2)*((2*P*x[i])-x[i]**2)
            dyc_dx[i] = ((2*M)/(P**2)*(P-x[i]))
        elif x[i] >= P and x[i] <= 1:
            yc[i] = (M/(1-P)**2)*(1-(2*P)+(2*P*x[i])-(x[i]**2))
            dyc_dx[i] = ((2*M)/((1-P)**2)*(P-x[i]))
        theta[i] = np.arctan(dyc_dx[i])
    
    yt = 5*T*((a0*x**0.5)+(a1*x) + (a2*x**2) + (a3*x**3) + (a4*x**4))
    xu = x - yt*np.sin(theta.T)
    yu = yc.T + yt*np.cos(theta.T)
    
    
    list_xu = xu.tolist()
    list_xu = list_xu[0]
    list_xu.reverse()
    list_yu = yu.tolist()
    list_yu = list_yu[0]
    list_yu.reverse()
    xu = np.array(list_xu)
    yu = np.array(list_yu)
    
    xl = x + yt*np.sin(theta.T)
    yl = yc.T - yt*np.cos(theta.T)
    
    fig = plt.figure(figsize=(8, 6))
    ax = fig.add_subplot(111)
    ax.plot(xu.T, yu.T, 'r')
    ax.plot(xl.T, yl.T, 'b')
    ax.set_aspect(1)
    
    return xu.T, yu.T, xl.T, yl.T


# Some basic geometrical data
halfSpan=1000.0
rootLength=100.0
tipLength=50.0
rootTwist=0.0
tipTwist=5.0


#connecting to windows com
CATIA = Dispatch('CATIA.Application')
CATIA.Visible = True

#create  an empty part
part_document = CATIA.Documents.Add("Part")
part1 = part_document.Part

#Shape factory provides generating of shapes
ShFactory = part1.HybridShapeFactory

#starting new body (geometrical set) in part1
bodies = part1.HybridBodies

#adding new body to part1
body1 = bodies.Add()
body1.Name = "Wireframe"

bodies2 = body1.hybridBodies  # Starting new geometrical set in Wireframe
body2 = bodies2.Add()  # Adding new body to Wireframe
body2.Name = "RootSection"  # Naming new body as "RootSection"

body3 = bodies2.Add()
body3.Name = "TipSection"

body4 = bodies.Add()  # Adding new body in part1
body4.Name = "Surfaces"  # Naming new body as "Surfaces"


# Creating new point [0,0,0] in Wireframe
point0 = ShFactory.AddNewPointCoord(0.000000, 0.000000, 0.000000)
body1.AppendHybridShape(point0)
# part1 should be updated after every new object
part1.Update() 

#Creatinging Z­direction for translating wing sections
wingAxis1= ShFactory.AddNewDirectionByCoord(0.000000, 0.000000, 1.000000)

#Creating twist point, sections will be twisted around this point
twistPoint1=ShFactory.AddNewPointCoord(25.0,0.0,0.0)
twistRef1= part1.CreateReferenceFromObject(twistPoint1)

#Creating Z­direction for translating wing sections
twistDir1 = ShFactory.AddNewDirectionByCoord(0.000000, 0.000000, 1.000000)
#Creating [POINT­DIRECTION] axis for twisting wing sections 
twistAxis1 = ShFactory.AddNewLinePtDir(twistRef1, twistDir1, 0.000000, 20.000000, False)

spline1 = ShFactory.AddNewSpline()
spline1.SetSplineType(0)
spline1.SetClosing(0)

xu, yu, xl, yl = NACA_Profil()

# Filling the spline with points
for i in range(0, len(xu)):
    
    point = ShFactory.AddNewPointCoord(xu[i]*100 , yu[i]*100, 50)  # coordinates are 2D, Z=0.0
    spline1.AddPoint(point)
    part1.Update()# new point to spline is added
for j in range(1, len(xl)):
    point2 = ShFactory.AddNewPointCoord(xl[j]*100 , yl[j]*100, 50)
    
    spline1.AddPoint(point2)
    part1.Update()
ShFactory.GSMVisibility(spline1,1)  # hide the spline
    
ref1 = part1.CreateReferenceFromObject(spline1)    
ref2 = part1.CreateReferenceFromObject(twistPoint1)
scaling1 = ShFactory.AddNewHybridScaling(ref1,ref2, rootLength/100.0) 
scaling1.VolumeResult = False
body2.AppendHybridShape(scaling1)
ShFactory.GSMVisibility(scaling1,1)

#Rotate [AXIS] the root section
rotate1= ShFactory.AddNewEmptyRotate()
ref1= part1.CreateReferenceFromObject(scaling1)
ref2 = part1.CreateReferenceFromObject(twistAxis1)
rotate1.ElemToRotate = ref1
rotate1.VolumeResult = False
rotate1.RotationType = 0
rotate1.Axis = twistAxis1
rotate1.AngleValue = rootTwist
body2.AppendHybridShape(rotate1)
ShFactory.GSMVisibility(rotate1,1)

#Translate [DIRECTION ­ DISTANCE] the root section
# is actually not necessary here
translate1 = ShFactory.AddNewEmptyTranslate()
ref1= part1.CreateReferenceFromObject(rotate1)
translate1.ElemToTranslate = rotate1
translate1.VectorType = 0
translate1.Direction = wingAxis1
translate1.DistanceValue = 0.00
translate1.VolumeResult = False
translate1.Name = "rootShape"
# Naming result "rootShape" IMPORTANT!!!
body2.AppendHybridShape(translate1)
part1.Update()

#%%-----------------------------------
#Creatinging Z­direction for translating wing sections
wingAxis2= ShFactory.AddNewDirectionByCoord(0.000000, 0.000000, 1.000000)

#Creating twist point, sections will be twisted around this point
twistPoint2=ShFactory.AddNewPointCoord(12.5,0.0,0.0)
twistRef2= part1.CreateReferenceFromObject(twistPoint2)

#Creating Z­direction for translating wing sections
twistDir2 = ShFactory.AddNewDirectionByCoord(0.000000, 0.000000, 1.000000)
#Creating [POINT­DIRECTION] axis for twisting wing sections 
twistAxis2 = ShFactory.AddNewLinePtDir(twistRef2, twistDir2, 0.000000, 20.000000, False)


spline2 = ShFactory.AddNewSpline()
spline2.SetSplineType(0)
spline2.SetClosing(0)


for i in range(0, len(xu)):
    
    point = ShFactory.AddNewPointCoord(xu[i]*50 , yu[i]*50, -50)  # coordinates are 2D, Z=0.0
    spline2.AddPoint(point)
    part1.Update()# new point to spline is added
for j in range(1, len(xl)):
    point2 = ShFactory.AddNewPointCoord(xl[j]*50 , yl[j]*50, -50)
    
    spline2.AddPoint(point2)
    part1.Update()
ShFactory.GSMVisibility(spline2,1)  # hide the spline
    
ref3 = part1.CreateReferenceFromObject(spline2)    
ref4 = part1.CreateReferenceFromObject(twistPoint2)
scaling2 = ShFactory.AddNewHybridScaling(ref3,ref4, tipLength/50.0) 
scaling2.VolumeResult = False
body3.AppendHybridShape(scaling2)
ShFactory.GSMVisibility(scaling2,1)

#Rotate [AXIS] the root section
rotate2= ShFactory.AddNewEmptyRotate()
ref3= part1.CreateReferenceFromObject(scaling2)
ref4 = part1.CreateReferenceFromObject(twistAxis2)
rotate2.ElemToRotate = ref3
rotate2.VolumeResult = False
rotate2.RotationType = 0
rotate2.Axis = twistAxis2
rotate2.AngleValue = rootTwist
body3.AppendHybridShape(rotate2)
ShFactory.GSMVisibility(rotate2,1)

#Translate [DIRECTION ­ DISTANCE] the root section
# is actually not necessary here
translate2 = ShFactory.AddNewEmptyTranslate()
ref3= part1.CreateReferenceFromObject(rotate2)
translate2.ElemToTranslate = rotate2
translate2.VectorType = 0
translate2.Direction = wingAxis1
translate2.DistanceValue = 0.00
translate2.VolumeResult = False
translate2.Name = "tipShape"
# Naming result "rootShape" IMPORTANT!!!
body3.AppendHybridShape(translate2)
part1.Update()

#Create new loft ­ MULTISECTION SURFACE
loft1 = ShFactory.AddNewLoft()
loft1.SectionCoupling = 1
loft1.Relimitation = 1
loft1.CanonicalDetection = 2
#Adding root section to the loft
shapes1 = body2.HybridShapes
# getting item from pool!!
result1 = shapes1.Item("rootShape")
ref1 = part1.CreateReferenceFromObject(result1)
ref2 = None
loft1.AddSectionToLoft(ref1, 1, ref2)
#Adding tip section to the loft
shapes2 = body3.HybridShapes
# getting item from pool!!
result2 = shapes2.Item("tipShape")
ref1 = part1.CreateReferenceFromObject(result2)
ref2 = None
loft1.AddSectionToLoft(ref1, 1, ref2)
loft1.Name = "masterSurface"
#Adding loft to Surfaces geometrical set
body4.AppendHybridShape(loft1)


part1.Update() 
    
    