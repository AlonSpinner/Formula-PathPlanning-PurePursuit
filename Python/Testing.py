# -*- coding: utf-8 -*-
"""
Created on Wed Apr 24 18:21:25 2019

@author: Alon
"""

import OrderConesDT
import matplotlib.pyplot as plt
import PurePursuit
import matplotlib.animation as animation

Json=OrderConesDT.JsonExpi('Json_18.txt')
Fig=plt.figure(figsize=(18,16))
Ax=Fig.add_subplot(111)
Ax.axis('equal')
MovieFreq=20.0 #Hz

StopFlag=True
for Itr in range(50,Json.FrameAmount):
    [Cones,CarCG,CarDir]=Json.GetFrame(Itr)

    MapTrack = OrderConesDT.MapTrack(Cones, CarCG, CarDir,SphereR=0.8,CarLength=0.1)  # Build class
    MapTrack.OrderCones(MaxItrAmnt=20, CostThreshold=-0.2, ColorCostWeight=0.4, \
               RRatioThreshold=10)
    MidPoints=MapTrack.FindMidPoints()
    MapTrack.PlotMap(Ax) #Map of cones   and all related

    PP=PurePursuit.PPAckerman(CarCG,CarDir,MidPoints,Lb=0.1,SymSteeringAngleBounds=12*3.1415/180,SphereR=0.8)
    PP.PlotMap(Ax) #plot PurePersuit related

    #last movie stuff
    Ax.set_title('frame number %g, steering angle %g' %(Itr,round((180/3.1415)*PP.SteeringAngle)),fontsize=24)

    plt.draw()
    plt.pause(1.0/MovieFreq)
    #raw_input("Press Enter to continue...") #activate for one by one
    plt.cla()