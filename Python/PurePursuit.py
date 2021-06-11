import numpy as np
import matplotlib.pyplot as plt

'''
Book:
Steering Control of an Autonomous Ground Vehicle with Application to the DARPA
Urban Challenge
By
Stefan F. Campbell
'''

class PPAckerman:
    def __init__(self,CarCG,CarDir,WayPoints,Lb=0.1,SymSteeringAngleBounds=12*3.1415/180,SphereR=1):
        self.CarCG=CarCG #[x,y]
        self.CarDir=CarDir #[x,y] - normalized
        self.WayPoints=WayPoints #mx2 [x,y] matrix of waypoints
        self.SymSteeringAngleBounds=SymSteeringAngleBounds #radians
        self.SphereR=SphereR #10 is an expiramental number
        self.Lb=Lb #distance between front and back wheels
        self.Goal=self.FindGoal() #[x,y] waypoint for which to aim. Initalized on startup
        self.SteeringAngle=self.FindSteeringAngle() #in radians. Initalized on startup

    def FindGoal(self):
        '''
        #Ackerman steering angle from page 87 (chapter 4.2, equation 4,2)
        :return: finds goal points from waypoints
        '''
        WayPoints = self.WayPoints

        #Remove WayPoints that  beyond CircumRadius of the first three points (around CarCG)
        if WayPoints.shape[0]>2:
            CircumR=np.amin([CircumRadius(WayPoints[:3,:]),self.SphereR])
            WayPoints=SphereFilter(WayPoints,self.CarCG,2*CircumR)[0]

        if WayPoints.size==0:
                Goal=self.CarCG+self.CarDir*2*self.Lb
        else:
            Goal=WayPoints[-1, :]  # pick last waypoint
        self.Goal=Goal
        return Goal

    def FindSteeringAngle(self):
        '''
        #Ackerman steering angle from page 87 (chapter 4.2, equation 4,2)
        :return: returns angles in radians
        '''
        La=np.linalg.norm(self.Goal-self.CarCG) #Look ahead distance
        n=(self.Goal-self.CarCG)/La
        eta=np.arctan2(np.cross(self.CarDir,n),np.dot(self.CarDir,n)) #angle of vector n relative to CarDir in normal convention
        SteeringAngle=np.arctan2(2*self.Lb*np.sin(eta),La) #Bicycle model

        if SteeringAngle>self.SymSteeringAngleBounds:
            SteeringAngle=self.SymSteeringAngleBounds
        if SteeringAngle<-self.SymSteeringAngleBounds:
            SteeringAngle=-self.SymSteeringAngleBounds

        self.SteeringAngle=SteeringAngle
        return SteeringAngle

    def PlotMap(self,Ax=plt.axes):
        #Ax.scatter(self.WayPoints[:, 0], self.WayPoints[:, 1], color=[0.5, 0, 0], s=70, # scatter waypoints
        #          edgecolors=[0.5, 0, 0])  # scatter waypoints
        #Ax.plot(self.WayPoints[:, 0], self.WayPoints[:, 1], color=[0.5, 0, 0], lw=2)  # plot waypoints interpolation

        Ax.scatter(self.Goal[0],self.Goal[1], color=[0.5, 0, 0], s=200, #scatter goal
                   facecolors='none' ,edgecolors=[0.5, 0, 0])

        #Ax.scatter(self.CarCG[0],self.CarCG[1], color=[0.5, 0, 0.5])  # plot CarCG
        #Ax.quiver(self.CarCG[0],self.CarCG[1], self.CarDir[0],self.CarDir[1])  # plot car direction

        v=RotateVector(self.CarDir,self.SteeringAngle)
        Ax.quiver(self.CarCG[0], self.CarCG[1],v[0],v[1], color=[1,0.8,0])  # plot wheels direction

        Ax.plot([self.CarCG[0],self.Goal[0]],[self.CarCG[1],self.Goal[1]], lw=0.5, color=[0,0.5,0], ls='dashed')

        Ax.grid(color='k', linestyle='-', linewidth=0.2)  # add grid
def SphereFilter(Cones, CarCG, R):
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
    ConesR = Cones[:, :2] - CarCG  # vectors of cones relative to car
    SqDistance = np.diag(np.matmul(ConesR, np.transpose(ConesR)))
    FInd = SqDistance < R ** 2
    FilteredCones = Cones[FInd, :]
    return FilteredCones, FInd
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


