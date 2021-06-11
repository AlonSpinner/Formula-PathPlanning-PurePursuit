# -*- coding: utf-8 -*-
import numpy as np
import matplotlib.pyplot as plt
from scipy.spatial import Delaunay
import json

class MapTrack:
    def __init__(self,Cones,CarCG,CarVel\
            ,SphereR=100, CarLength=0.1):
        #initial geomteric mapping
        self.Cones=Cones #matrix mx3 [x,y,color]. color = 98(blue)/121(yellow)
        self.CarCG=CarCG #[x,y]
        self.CarDir=CarVel/np.linalg.norm(CarVel) #[x,y] of car direction

        #Geomtetric parameters
        self.SphereR=SphereR #10 is an expiramental number
        self.CarLength=CarLength#0.1 is an expiramental number

        #initalize output of algorithm
        self.Cones4Calc=self.FindCones4Calc() #cones from which track is calculated.
        self.OrderedBlueCones=np.array([]) #Calculated in OrderCones
        self.OrderedYellowCones=np.array([]) #Calculated in OrderCones
        self.MidPoints=np.array([]) #Calculated in FindMidPoints
    def FindCones4Calc(self,Angle4FakeCones=np.pi/3):
        #Preprocessing
        Cones=SphereFilter(self.Cones,self.CarCG,self.SphereR)[0] #filtering for radius around car

        #Fix new points - fake cones - if no blue/yellow cones are found
        BlueCones=Cones[(Cones[:,2]==98),:];  FrontBlueCones=BackFilter(BlueCones, self.CarCG, self.CarDir)[0]
        YellowCones = Cones[(Cones[:, 2] == 121), :]; FrontYellowCones = BackFilter(YellowCones, self.CarCG, self.CarDir)[0]
        BaddFlag=FrontBlueCones.size==0
        YaddFlag=FrontYellowCones.size==0
        if BaddFlag:
            BNew=self.CarCG+self.CarLength*RotateVector(self.CarDir,+Angle4FakeCones) #add cone to the left
            Cones=np.vstack([Cones,np.hstack([BNew,98])]) #add fake yellow cone @ the buttom of matrix
        if YaddFlag: #no yellow cones exist
            YNew=self.CarCG+self.CarLength*RotateVector(self.CarDir,-Angle4FakeCones) #add cone to the right
            Cones=np.vstack([Cones,np.hstack([YNew,121])]) #add fake yellow cone @ the buttom of matrix

        self.Cones4Calc=Cones
        return Cones
    def OrderCones(self,MaxItrAmnt=10,CostThreshold=-0.2,ColorCostWeight=0.6, \
                   RRatioThreshold=20):
        '''
        Input:
        self

        Output:
        to class:
        BCones - blue cones in order for interpolation [x,y]
        Ycones - yellow cones in order for interpolation [x,y]
        '''
        if self.Cones4Calc.shape[0] < 3:
            return 0

        #initalize
        Yind=np.full([MaxItrAmnt,1],np.nan)
        Bind=np.full([MaxItrAmnt,1],np.nan)
        Cones=self.Cones4Calc

        #Triangulate
        DT=Delaunay(Cones[:,:2]) #Triangulate only /w Cones
        ID=DT.find_simplex(self.CarCG) #attempt to find triangle which contains CarCG
        if ID==-1: #if CarCG isoutside of convex hull
            Cones=np.vstack([Cones,np.hstack([self.CarCG,0])]) #add Car to Cones as a fake cone
            DT=Delaunay(Cones[:,:2])#Triangulate /w Cones+CarCG
            ID=DT.find_simplex(self.CarCG+0.5*self.CarLength*self.CarDir) #find first triangle to work with
        if ID == -1: return  # no cones in sight. couldnt build a track

        #Important note: DT.find_simplex returns -1 if point is outside triangulation.
        #But we use setdiff1d in both TriOne and FindNextTriangle so NewID can be of empty.

        NewID,Newu,NewCrossEdge=TriOne(DT,Cones[:,2],ID,self.CarDir,ColorCostWeight)
        if NewID.size==0: return #crossed into no-man's land

        #insert first cones into lists (Bind/Yind)
        if Cones[NewCrossEdge[0],2]==ord('y'):
            Yind[0]=NewCrossEdge[0]
        else:
            Bind[0]=NewCrossEdge[0]
        if Cones[NewCrossEdge[1],2]==ord('y'):
            Yind[0]=NewCrossEdge[1]
        else:
            Bind[0]=NewCrossEdge[1]
        NewV=np.setdiff1d(DT.simplices[NewID,:],NewCrossEdge);
        if Cones[NewV,2]==ord('y'):
            Yind[1]=NewV
        else:
            Bind[1]=NewV

        for Itr in range(2,MaxItrAmnt):
            #find next triangle
            NewID,Newu,NewCrossEdge,Cost,RRatio=FindNextTriangle(DT,Cones[:,2],NewID,Newu,NewCrossEdge,ColorCostWeight)
            #check conditions, if not good enough - break
            if NewID.size==0 or CostThreshold<Cost or RRatioThreshold<RRatio:
                break
            #add next cone to Bind/Yind
            NewV=np.setdiff1d(DT.simplices[NewID,:],NewCrossEdge)
            if Cones[NewV,2]==ord('y'):
                Yind[Itr]=NewV
            else:
                Bind[Itr]=NewV

        #create output
        Bind=Bind[~np.isnan(Bind)].astype(int); Yind=Yind[~np.isnan(Yind)].astype(int) #delete rows with NaNs
        BCones=Cones[Bind,:2]; YCones=Cones[Yind,:2];

        #Output to class
        self.OrderedBlueCones=BCones
        self.OrderedYellowCones=YCones
    def FindMidPoints(self):
        '''
        Input:
            self
        Output:
            MiddlePoints - ordered points of middle of track [x,y]
        '''
        #decide on inner and outer cones by amount (infantile)
        Bamnt=self.OrderedBlueCones.shape[0]; Yamnt=self.OrderedYellowCones.shape[0]
        Mamnt=np.maximum(Bamnt,Yamnt) #amount of middle points
        Ind=np.argmax([Bamnt,Yamnt]) #index of maximum between [Bamnt,Yamnt]
        if Ind==0:
            Morexy=self.OrderedBlueCones; Lessxy=self.OrderedYellowCones
        else:
            Morexy=self.OrderedYellowCones; Lessxy=self.OrderedBlueCones

        #Build the middle points
        MidPoints=np.zeros([Mamnt,2])
        for k in range(0,Mamnt):
            MoreCone=Morexy[k,:]
            Ind=np.argmin(np.diag(np.matmul(Lessxy-MoreCone,np.transpose(Lessxy-MoreCone)))) #find closest InCone to k-th OutCone
            SmallCone=Lessxy[Ind,:]
            MidPoints[k,:]=(SmallCone+MoreCone)/2

        #Remove Points behind the vehicle
        MidPoints=BackFilter(MidPoints,self.CarCG,self.CarDir)[0]

        self.MidPoints=MidPoints
        return MidPoints
    def PlotMap(self,Ax=plt.axes):
        '''
        :param Ax: input axes to draw upon. if none exist function will create.
        :return: figure with scattered cones, interpolated cones, midpoints, carCG, car direction and lookahead distance.
        '''

        #Initialize
        BlueCones = self.Cones[self.Cones[:, 2] == 98] #extract array of blue cones
        YellowCones = self.Cones[self.Cones[:, 2] == 121] #extract array of yellow cones
        MidPoints = self.MidPoints #extract midpoints

        #Plot it all up
        Ax.scatter(BlueCones[:, 0], BlueCones[:, 1], color=[0, 0, 1], s=70,
                    edgecolors=[0, 0, 1])  # scatter blue cones
        Ax.scatter(YellowCones[:, 0], YellowCones[:, 1], color=[0.9, 0.9, 0], s=70,
                    edgecolors=[0.9, 0.9, 0])  # scatter yellow cones
        Ax.scatter(self.Cones4Calc[:, 0], self.Cones4Calc[:, 1], color=[0, 0, 0], s=100,
                   facecolors='none')  # plot Cones4Calc
        if self.OrderedBlueCones.size>0:
            plt.plot(self.OrderedBlueCones[:, 0], self.OrderedBlueCones[:, 1], color=[0, 0, 1],
                 lw=2)  # plot yellow cones interpolation
        if self.OrderedYellowCones.size>0:
            plt.plot(self.OrderedYellowCones[:, 0], self.OrderedYellowCones[:, 1], color=[0.7, 0.7,0],
                 lw=2)  # plot blue cones interpolation
        if self.MidPoints.size > 0:
            Ax.scatter(MidPoints[:, 0], MidPoints[:, 1], color=[0.5, 0, 0], s=70,
                        edgecolors=[0.5, 0, 0])  # scatter midpoints
            Ax.plot(MidPoints[:, 0], MidPoints[:, 1], color=[0.5, 0, 0], lw=2)  # plot middle points interpolation
        Ax.scatter(self.CarCG[0],self.CarCG[1], color=[0.5, 0, 0.5])  # plot CarCG
        Ax.quiver(self.CarCG[0],self.CarCG[1], self.CarDir[0],self.CarDir[1])  # plot car direction

        t = np.linspace(0, 2 * np.pi, 100)  # for circle plotting
        Ax.plot(self.CarCG[0]+self.SphereR*np.cos(t),self.CarCG[1]+self.SphereR*np.sin(t),ls='--',lw=2,color=[0,0.5,0]) #plot circle
        Ax.grid(color='k', linestyle='-', linewidth=0.2)  # add grid
#Associated functions with class ConesDT
def TriOne(DT,ConesColors,ID,CarDir,ColorCostWeight):
    '''
    Input:
    DT - DelaunayTriangulation containing DT.points and DT.simplices
    ConesColors - colors of cones (98 for blue, 121 for yellow) with same
    indexing as DT.Points
    ID - number of triangle in DT.simplices (row index)
    CarDir - car direction [x,y] normalized
    ColorCostWeight - weight for costfunction in triangles. Applied on first term - effect cone colors.

    Output:
    NewID - number of new triangle in DT.ConnectivityList (row index)
    Newu - direction of enterance to new triangle (NewID)
    NewCrossEdge - [V1,V2] of edge that we crossed to get from ID->NewID
    V1 and V2 refer to vertex indcies (row) of DT.points
    '''
    V=DT.simplices[ID,:] #find vertex indcies of ID
    P=np.hstack([(DT.points[V,:]),ConesColors[V].reshape(3,1)]) #Triangle cones
    IC=InCenter(P[:,:2])[0] #find incenter of ID
    EdgesInd=np.array([[0,1], #Edge1 #build indcies in V mapping to EdgesV by V[EdgesInd]
                       [0,2], #Edge2
                       [1,2]])#Edge3
    J1=EdgeCost(IC,P[EdgesInd[0,:],:],CarDir,ColorCostWeight)
    J2=EdgeCost(IC,P[EdgesInd[1,:],:],CarDir,ColorCostWeight)
    J3=EdgeCost(IC,P[EdgesInd[2,:],:],CarDir,ColorCostWeight)
    J=np.hstack([J1,J2,J3]); JminInd=np.argmin(J)
    NewCrossEdge=V[EdgesInd[JminInd,:]] #obtain Edge Vertices to cross (point indcies [V1,V2])
    Newu=OutFacingNormal(P[EdgesInd[JminInd,:],:2],IC) #calculate normal to cross edge
    Attachments=np.where(np.any(DT.simplices==NewCrossEdge[0],axis=1) &\
                         np.any(DT.simplices==NewCrossEdge[1],axis=1)) #find IDs that are connected to edge
    NewID=np.squeeze(np.setdiff1d(Attachments,ID)) #NewID - New Triangle
    return NewID,Newu,NewCrossEdge
def FindNextTriangle(DT,ConesColors,ID,Dir,CrossEdge,ColorCostWeight):
    '''
    Input:
    DT - DelaunayTriangulation containing DT.points and DT.simplices
    ConesColors - colors of cones (98 for blue, 121 for yellow) with same
    indexing as DT.Points
    ID - number of triangle in DT.ConnectivityList (row index)
    Dir - direction of enterance to triangle ID
    CrossEdge - [V1,V2] of edge that was crossed to enter triangle ID
    V1 and V2 refer to vertex indcies (row) in DT.Points
    ColorCostWeight - weight for costfunction in triangles. Applied on first term - effect cone colors.

    Output:
    NewID - number of new triangle in DT.ConnectivityList (row index)
    NewDir - direction of enterance to new triangle (NewID)
    NewCrossEdge - [V1,V2] of edge that we crossed to get from ID->NewID
    V1 and V2 refer to vertex indcies (row) in DT.Points
    MinCost - price in cost function with which algorithm decided to go
    RRatio - ratio between Circumcenter/InCenter radii
    '''
    V=DT.simplices[ID,:] #find vertex indcies of ID
    P=np.hstack([DT.points[V,:],ConesColors[V].reshape(3,1)]) #Triangle cones
    IC=InCenter(P[:,:2])[0] #find incenter of ID
    EdgesInd=np.array([[0,1], #Edge1 #build indcies in V mapping to EdgesV by V[EdgesInd]
                       [0,2], #Edge2
                       [1,2]])#Edge3
    EdgesV=V[EdgesInd]
    EdgesInd=EdgesInd[~(np.any(EdgesV==CrossEdge[0],axis=1) &\
                        np.any(EdgesV==CrossEdge[1],axis=1)),:] #find the two edges to calculate for (not going backwards).
    J1=EdgeCost(IC,P[EdgesInd[0,:],:],Dir,ColorCostWeight) #calculate costs
    J2=EdgeCost(IC,P[EdgesInd[1,:],:],Dir,ColorCostWeight)
    J=np.hstack([J1,J2]); MinCost=np.min(J); JminInd=np.argmin(J) #find Edge to cross by row index in Edges
    NewCrossEdge=V[EdgesInd[JminInd,:]] #obtain Edge Vertices to cross (point indcies [V1,V2])
    NewDir=OutFacingNormal(P[EdgesInd[JminInd,:],:2],IC) #calculate normal to cross edge
    Attachments=np.where(np.any(DT.simplices==NewCrossEdge[0],axis=1) &\
                         np.any(DT.simplices==NewCrossEdge[1],axis=1)) #find IDs that are connected to edge
    NewID=np.squeeze(np.setdiff1d(Attachments,ID)) #NewID - New Triangle

    if NewID.size==0:
        RRatio=0 #filler
        return NewID,NewDir,NewCrossEdge,MinCost,RRatio

    NewV=DT.simplices[NewID,:] #find vertex indcies of NewID (=NewTriangle)
#    NewV=NewV.reshape(3)
    NewP=DT.points[NewV,:] #find points correlating to new vertex indcies
    RRatio=CircumRadius(NewP)/InCenter(NewP)[1] #calculate RRatio of new triangle
    return NewID,NewDir,NewCrossEdge,MinCost,RRatio
def EdgeCost(TriInCenter,EdgeCones,u,ColorCostWeight):
    '''
    Input:
    InCenter - [x,y] in center of triangle
    EdgeCones - [x,y,color] 2x3 of edge vertcies
    u - [x,y] normalized direction of enterance to triangle
    ColorCostWeight - weight for different color cones cost. should range [0,1]

    Output:
    J - cost of passing through edge. |J|<2

    Example:
    t1=np.radians([-30,90,210])
    ConesColors=np.array([98,121,98])
    TriCones=np.transpose(np.vstack([(np.cos(t1)),(np.sin(t1)),ConesColors])) #equilateral triangle
    EdgeCones=TriCones[:2,:] #edge cones have different colors
    InCenter=InCenter(TriCones[:,:2])[0] #[0] after function call - returns the  first value
    u=OutFacingNormal(EdgeCones[:,:2],InCenter) #enterance direction - same direction as out facing normal of edge
    CostColorWeight=1
    J=EdgeCost(InCenter,EdgeCones,u,alpha)
    print(J) #expect to be -2
    '''
    J=0 #initalize
    if EdgeCones[0,2]==EdgeCones[1,2]: #if the same color
        J=J+1
    else: #not the same color
        J=J-1
    J=ColorCostWeight*J-np.dot(u,OutFacingNormal(EdgeCones[:,:2],TriInCenter)) #+bad points for difference in direction
    return J
def BackFilter(Cones,CarCG,CarDir):
        '''
        Input:
        Cones - mx3 cones matrix containing [x,y,color]
        CarCG - [x,y] of car position
        CarDir - [x,y] normalized vector of car direction

        Output:
        FilteredCones - same format as cones. only cones that are infront of the vehicle
        FInd - bool array FilteredCones=Cones(FInd,:)
        '''
        ConesR=Cones[:,:2]-CarCG #vectors of cones relative to car
        DotProducts = np.matmul(ConesR, np.transpose(CarDir))  # dot product with normal (car direction)
        FInd = DotProducts > 0  # find MidPoints infront of the vehicle (index)
        FilteredCones=Cones[FInd,:]
        return FilteredCones,FInd
def SphereFilter(Cones,CarCG,R):
        '''
       Input:
        Cones - mx3 cones matrix containing [x,y,color]
        CarCG - [x,y] of car position
        R - radius of hemisphere

        Output:
        FilteredCones - same format as cones
        only cones within sphere with radius R starting at car
        FInd - bool array FilteredCones=Cones(FInd,:)

        Example:
        R=np.sqrt(2)
        Cones=np.array([[-1,0.5,98],
                        [1,1,121],
                        [1,-0.6,98]])
        CarCG=np.array([0,0.5])
        CarDir=np.array([0,1])
        SphereFilter(Cones,CarCG,R)
        PassCones, FInd=SphereFilter(Cones,CarCG,R) #run function on data
        FailCones=Cones[~FInd,:] #obtain fail cones
        t=np.linspace(0,2*np.pi,100) #for circle plotting
        plt.scatter(PassCones[:,0],PassCones[:,1],color=[0,1,0]) #plot pass cones
        plt.scatter(FailCones[:,0],FailCones[:,1],color=[1,0,0]) #plot fail cones
        plt.scatter(CarCG[0],CarCG[1],color=[0.5,0,0.5]) #plot CarCG
        plt.plot(CarCG[0]+R*np.cos(t),CarCG[1]+R*np.sin(t),ls='--',lw=2,color=[0,0,0]) #plot circle
        plt.grid(color='k',linestyle='-',linewidth=0.2) #add grid
        '''
        ConesR=Cones[:,:2]-CarCG #vectors of cones relative to car
        SqDistance=np.diag(np.matmul(ConesR,np.transpose(ConesR)))
        FInd=SqDistance < R**2
        FilteredCones=Cones[FInd,:]
        return FilteredCones,FInd

#Some basic geometrey functions we need
def ProjPnt2Line(P,q):
    '''
    Input:
    q - [x,y] point to be project onto line P
    P - 2x2 [[x1,y1], points representing a line
            [x2,y2]]

    Output:
    projq  [x,y] np array representing point resulting in projecting q on P

    Example:
    P=np.array([[0,0],
              [1,0]])
    q=np.array([0.5,4])
    projq=ProjPnt2Line(P,q)
    plt.scatter(P[:,0],P[:,1],color=[1,0,0])
    plt.scatter(q[0],q[1],color=[0,0,0])
    plt.scatter(projq[0],projq[1],color=[0,0,1])
    '''
    p1=P[0,:]; p2=P[1,:]
    t=(p2-p1)/np.linalg.norm(p2-p1) #unit vector in p1-p2
    projq=np.dot((q-p1),t)*t+p1; #find project point
    return projq
def InCenter(P):
    '''
    Input:
    P - 3x2 [[x1,y1], points of a triangle
             [x2,y2],
             [x3,y3]]

    Output:
    TriInCenter - [x,y] coordinates
    InR - radius of InCircle

    From:
    https://www.mathopenref.com/coordincenter.html
    https://keisan.casio.com/

    Example:
    t1=np.radians([-30,90,210])
    P=np.transpose(np.vstack([(np.cos(t1)),(np.sin(t1))])) #equilateral triangle
    InCenter, InR=InCenter(P)
    plt.plot(np.hstack([P[:,0],P[0,0]]),np.hstack([P[:,1],P[0,1]]),color=[0.5,0,0]) #plot triangle edges
    plt.scatter(P[:,0],P[:,1],color=[1,0,0]) #plot triangle points
    plt.scatter(InCenter[0],InCenter[1],color=[0,0,1]) #plot InCenter
    t2=np.linspace(0,2*np.pi,100) #for circle plotting
    plt.plot(InCenter[0]+InR*np.cos(t2),InCenter[1]+InR*np.sin(t2),ls='--',lw=2,color=[0,0,0]) #plot InCircle
    q=InCenter+InR*(P[0,:]-InCenter)/np.linalg.norm((P[0,:]-InCenter)) #find projection of point 1 on InCircle
    plt.plot([InCenter[0],q[0]],[InCenter[1],q[1]],lw=2,color=[0,0.7,0]) #plot in radius in the direction of point 1
    plt.grid(color='k',linestyle='-',linewidth=0.2) #add grid
    '''
    A=P[0,:]; B=P[1,:]; C=P[2,:]
    a=np.linalg.norm(B-C); b=np.linalg.norm(A-C); c=np.linalg.norm(A-B);
    p=a+b+c
    TriInCenter=(a*A+b*B+c*C)/p
    s=p/2
    InR=np.sqrt(s*(s-a)*(s-b)*(s-c))/s
    return TriInCenter,InR
def CircumRadius(P):
    '''
    Input:
    P - 3x2 [[x1,y1], points of a triangle
             [x2,y2],
             [x3,y3]]

    Output:
    CircumR=radius of CircumCircle

    From:
    http://mathworld.wolfram.com/Circumradius.html

    Example:
    t1=np.radians([-30,90,210])
    P=np.transpose(np.vstack([(np.cos(t1)),(np.sin(t1))])) #equilateral triangle
    CircumCenter=[0,0] #symmetry
    CircumR=CircumRadius(P)
    plt.plot(np.hstack([P[:,0],P[0,0]]),np.hstack([P[:,1],P[0,1]]),color=[0.5,0,0]) #plot triangle edges
    plt.scatter(P[:,0],P[:,1],color=[1,0,0]) #plot triangle points
    plt.scatter(CircumCenter[0],CircumCenter[1],color=[0,0,1]) #plot InCenter
    t2=np.linspace(0,2*np.pi,100) #for circle plotting
    plt.plot(CircumCenter[0]+CircumR*np.cos(t2),CircumCenter[1]+CircumR*np.sin(t2),ls='--',lw=2,color=[0,0,0]) #plot InCircle
    plt.grid(color='k',linestyle='-',linewidth=0.2) #add grid
    '''
    A=P[0,:]; B=P[1,:]; C=P[2,:]
    a=np.linalg.norm(B-C); b=np.linalg.norm(A-C); c=np.linalg.norm(A-B);
    CircumR=a*b*c/np.sqrt((a+b+c)*(a+b-c)*(a+c-b)*(b+c-a))
    return CircumR
def OutFacingNormal(P,TriInCenter):
    '''
    Input:
    TriInCenter [x,y] of a triangle incenter
    P - 2x2 [[x1,y1], points on triangle of the same edge
             [x2,y2]]

    Output:
    n - unit vector [x,y] normal to the edge specified by P, facing out of the triangle

    Example:
    t1=np.radians([-30,90,210])
    Ptri=np.transpose(np.vstack([(np.cos(t1)),(np.sin(t1))])) #equilateral triangle
    Pedge=Ptri[:2,:]
    InCenter,InR=InCenter(Ptri)
    q=np.mean(Pedge,axis=0) #mean by column
    n=OutFacingNormal(Pedge,InCenter)
    plt.plot(np.hstack([Ptri[:,0],Ptri[0,0]]),np.hstack([Ptri[:,1],Ptri[0,1]]),color=[0.5,0,0]) #plot triangle
    plt.scatter(Pedge[:,0],Pedge[:,1],color=[1,0,0]) #plot edge points
    plt.scatter(InCenter[0],InCenter[1],color=[0,0,1]) #plot InCenter
    plt.quiver(q[0],q[1],n[0],n[1])
    plt.grid(color='k',linestyle='-',linewidth=0.2) #add grid
    '''
    q=ProjPnt2Line(P,TriInCenter);
    n=(q-TriInCenter)/np.linalg.norm(q-TriInCenter);
    return n
def RotateVector(n,Theta):
    '''
    Input:
        n - row vector [x,y]
        Theta - angle in radians to rotate by
    Output:
        v - row rotated vector [x,y]
    '''
    Q=np.array([[np.cos(Theta),-np.sin(Theta)],
               [np.sin(Theta),np.cos(Theta)]]) #Rotation matrix
    v=np.transpose(np.matmul(Q,np.transpose(n)))
    return v

#Importing Json files extracted from ROS to pything enviorment
class JsonExpi:
    '''
    Class for importing data coming out of ROS.
    Tom has configured printing to json format which is utilized by this class


    EXAMPLE:
    import JsonImport
    Data=JsonImport.JsonExpi('Json_03062019.txt')
    for FrameNum in range(1,Data.FrameAmount):
    [Cones,CarCG,CarDir]=Data.GetFrame(FrameNum)
    print(Cones); print(CarCG); print(CarDir)
    '''
    def __init__(self,FileName):
        self.FileName=FileName #string

        with open(FileName, 'r') as json_file:
            Parsed = [json.loads(line) for line in json_file]
        self.Parsed=Parsed #list of frames

        self.FrameAmount=len(Parsed) #amount of frames in list
    def GetFrame(self,FrameNum):
        Frame=self.Parsed[FrameNum-1][str(FrameNum)][0]

        Bx=np.array(Frame["blue_x"])
        By=np.array(Frame["blue_y"])
        Yx=np.array(Frame["yellow_x"])
        Yy=np.array(Frame["yellow_y"])
        CarDir=np.array(Frame["normal"])
        CarCG=np.array(Frame["pose"])

        BCones=np.transpose(np.vstack([Bx,By,98*np.ones_like(Bx)]))
        YCones = np.transpose(np.vstack([Yx, Yy, 121*np.ones_like(Yx)]))
        Cones=np.vstack([BCones,YCones])

        return Cones,CarCG,CarDir