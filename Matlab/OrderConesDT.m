function [BCones,YCones]=OrderConesDT(Cones,CarCG,CarVel)
%Input:
%Cones - double matrix 3xm [x,y,color]. color = 98(blue)/121(yellow)
%CarCG - double [x,y] of geometric center of car
%CarDir - double [x,y] of car velocity (direcational)

%Output:
%BCones - blue cones in order for interpolation [x,y]
%Ycones - yellow cones in order for interpolation [x,y]

%DelaunayTriangulation mathworks library
%https://www.mathworks.com/help/matlab/ref/delaunaytriangulation.html
%Some info about functions
%ID=pointLocation(DT,Q); %returns ConnectivityList row index of the triangle binding query point
%V=DT.ConnectivityList(ID,:); %returns the vertex indcies of the points making the triangle pointed by ID
%P=DT.Points(V,:); %returns the points of the triangle with vertices V
%[IC,r]=incenter(DT,ID); %center of circle whose bounded by the triangle
%[OC,R]=circumcenter(DT,NewID); %center of circle who binds the triangle

%running parameters
MaxItrAmnt=200; 
CostThreshold=-0.2; %if Cost>CostThreshold, exit
RRatioThreshold=10; %if Rbinding/Rbounded>RRatioThreshold, exit
SphereR=15000; %m
CarLength=5; %m

%initalize
Yind=nan(MaxItrAmnt,1);
Bind=nan(MaxItrAmnt,1);

%Preprocessing
CarDir=CarVel/norm(CarVel); %normalize velocity to unit direction
[~,FInd]=SphereFilter(Cones(:,[1,2]),CarCG,SphereR); %filtering for radius around car
Cones=Cones(FInd,:);

%Triangulate
DT=delaunayTriangulation(Cones(:,[1,2])); %Triangulate only /w Cones
ID=pointLocation(DT,CarCG); %attempt to find triangle which contains CarCG
if isnan(ID) %if CarCG isoutside of convex hull
    Cones=[Cones;[CarCG,0]]; %add Car to Cones as a fake cone
    DT=delaunayTriangulation(Cones(:,1:2)); %Triangulate /w Cones+CarCG
    ID=pointLocation(DT,CarCG+0.5*CarLength*CarDir); %find first triangle to work with  
end
[NewID,Newu,NewCrossEdge]=TriOne(DT,Cones(:,3),ID,CarDir);
if isempty(NewID), return, end %crossed into no-man's land

%insert first cones into lists (Bind/Yind)
if Cones(NewCrossEdge(1),3)=='y', Yind(1)=NewCrossEdge(1);
else, Bind(1)=NewCrossEdge(1); end
if Cones(NewCrossEdge(2),3)=='y', Yind(1)=NewCrossEdge(2);
else, Bind(1)=NewCrossEdge(2); end
NewV=setdiff(DT.ConnectivityList(NewID,:),NewCrossEdge);
if Cones(NewV,3)=='y', Yind(2)=NewV;
else, Bind(2)=NewV; end

for Itr=3:MaxItrAmnt    
    %find next triangle
    [NewID,Newu,NewCrossEdge,Cost,RRatio]=FindNextTriangle(DT,Cones(:,3),NewID,Newu,NewCrossEdge);
    %check conditions, if not good enough - break
    if isempty(NewID) || CostThreshold<Cost || RRatioThreshold<RRatio, break; end
    %add next cone to Bind/Yind
    NewV=setdiff(DT.ConnectivityList(NewID,:),NewCrossEdge);
    if Cones(NewV,3)=='y', Yind(Itr)=NewV;
    else, Bind(Itr)=NewV; end
end

%create output
Bind(any(isnan(Bind),2),:)=[]; Yind(any(isnan(Yind),2),:)=[]; %delete rows with NaNs
BCones=Cones(Bind,1:2); YCones=Cones(Yind,1:2);
end
function [NewID,Newu,NewCrossEdge]=TriOne(DT,ConesColors,ID,CarDir)
%Input:
%DT - DelaunayTriangulation containing DT.Points and DT.ConnectivityList
%ConesColors - colors of cones (98 for blue, 121 for yellow) with same
%indexing as DT.Points
%ID - number of triangle in DT.ConnectivityList (row index)
%CarDir - car direction [x,y] normalized

%Output:
%NewID - number of new triangle in DT.ConnectivityList (row index)
%Newu - direction of enterance to new triangle (NewID)
%NewCrossEdge - [V1,V2] of edge that we crossed to get from ID->NewID
%V1 and V2 refer to vertex indcies (row) in DT.Points

V=DT.ConnectivityList(ID,:); %find vertex indcies of ID
InCenter=incenter(DT,ID); %find incenter of ID

%build P=[x,y,color] of triangle points - "local triangle points [1,2,3]"
P=zeros(3,3);
P(:,1:2)=DT.Points(V,:); %obtain the triagnle points themselves
P(:,3)=ConesColors(V);

%find the two edges to calculate for (not going backwards).
Edges=[V(1),V(2);V(1),V(3);V(2),V(3)]; %build all edges
[~,Vind]=ismember(Edges,V); %maps V in matrix Edges to index numbers in V, and in consequence, in P
%V=[2,1,4]; Edges=[2,1;2,4]; ->Vind=[1,2;1,3];
%find EdgeCosts
J1=EdgeCost(InCenter,P(Vind(1,:),:),CarDir);
J2=EdgeCost(InCenter,P(Vind(2,:),:),CarDir);
J3=EdgeCost(InCenter,P(Vind(3,:),:),CarDir);
J=[J1,J2,J3];

%decide on edge to cross and calculate normal to it
[~,JInd]=min(J);
NewCrossEdge=Edges(JInd,:);
Newu=OutFacingNormal(P(Vind(JInd,:),[1,2]),InCenter);

%New ID of triangle
Attchments=edgeAttachments(DT,NewCrossEdge); %cell containing 2 triangles ID {current,next}
NewID=setdiff(Attchments{1},ID);
end
function [NewID,Newu,NewCrossEdge,MinCost,RRatio]=FindNextTriangle(DT,ConesColors,ID,u,CrossEdge)
%Input:
%DT - DelaunayTriangulation containing DT.Points and DT.ConnectivityList
%ConesColors - colors of cones (98 for blue, 121 for yellow) with same
%indexing as DT.Points
%ID - number of triangle in DT.ConnectivityList (row index)
%u - direction of enterance to triangle ID
%CrossEdge - [V1,V2] of edge that was crossed to enter triangle ID
%V1 and V2 refer to vertex indcies (row) in DT.Points

%Output:
%NewID - number of new triangle in DT.ConnectivityList (row index)
%Newu - direction of enterance to new triangle (NewID)
%NewCrossEdge - [V1,V2] of edge that we crossed to get from ID->NewID
%V1 and V2 refer to vertex indcies (row) in DT.Points
%MinCost - price in cost function with which algorithm decided to go
%RRatio - ratio between Circumcenter/InCenter radii

V=DT.ConnectivityList(ID,:); %find vertex indcies of ID
InCenter=incenter(DT,ID); %find incenter of ID

%build P=[x,y,color] of triangle points - "local triangle points [1,2,3]"
P=zeros(3,3);
P(:,1:2)=DT.Points(V,:); %obtain the triagnle points themselves
P(:,3)=ConesColors(V);

%find the two edges to calculate for (not going backwards).
Edges=[V(1),V(2);V(1),V(3);V(2),V(3)]; %build all edges
Edges=Edges(~ismember(Edges,[CrossEdge;fliplr(CrossEdge)],'rows'),:); %delete common edge with CrossEdge
[~,Vind]=ismember(Edges,V); %maps V in matrix Edges to index numbers in V, and in consequence, in P
%V=[2,1,4]; Edges=[2,1;2,4]; ->Vind=[1,2;1,3];
%find EdgeCosts
J1=EdgeCost(InCenter,P(Vind(1,:),:),u);
J2=EdgeCost(InCenter,P(Vind(2,:),:),u);
J=[J1,J2];

%decide on edge to cross and calculate normal to it
[MinCost,JInd]=min(J);
NewCrossEdge=Edges(JInd,:);
Newu=OutFacingNormal(P(Vind(JInd,:),[1,2]),InCenter);

%New ID of triangle
Attchments=edgeAttachments(DT,NewCrossEdge); %cell containing 2 triangles ID {current,next}
NewID=setdiff(Attchments{1},ID);

%If we crossed to "no man's land" - return
RRatio=[];
if isempty(NewID), return, end
%CircleRatio
[~,R]=circumcenter(DT,NewID);
[~,r]=incenter(DT,NewID);
RRatio=R/r;
end
function J=EdgeCost(InCenter,P,u)
%InCenter - in center of triangle [x,y]
%P=[x,y,color] of edges vertices 2x3
%u - direction from previous triangle to current triangle
J=0;
if abs((P(1,3)-P(2,3)))<eps %if they are the same color
    J=J+1; else, J=J-1; end
J=J-dot(u,OutFacingNormal(P(:,1:2),InCenter)); %+bad points for difference in direction
end
function n=OutFacingNormal(P,InCenter)
%function returns the out facing normal for an edge of a triangle [p1,p2]
%with InCenter.
%P=[x,y] of edges vertices 2x2
q=PorjPnt2Line(P(1,:),P(2,:),InCenter);
n=(q-InCenter)/norm(q-InCenter);
end
function projq=PorjPnt2Line(p1,p2,q)
%p1,p2,q,projq - [X,Y] points.
%p1,p2 - represent a line
t=(p2-p1)/norm(p2-p1);
projq=dot((q-p1),t)*t+p1;
end
function [FilteredCones,FInd]=SphereFilter(Cones,CarCG,R)
%Input:
%Cones - mx3 matrix of cones [x,y,color]
%CarCG - car center [x,y]
%CarDir - car direction [x,y] - normalized
%R - radius of hemisphere

%Output:
%FilteredCones - same format as cones
%only cones within sphere with radius R starting at car and openning in the Car's
%direction
%FInd - bool array FilteredCones=Cones(FInd,:);

ConesR=Cones(:,1:2)-CarCG;
SqDistance=diag(ConesR*ConesR');
FInd=(SqDistance<R^2);
FilteredCones=Cones(FInd,:);
end