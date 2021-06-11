function varargout = GenRaceway_limited(varargin)
% GENRACEWAY_LIMITED MATLAB code for GenRaceway_limited.fig
%      GENRACEWAY_LIMITED, by itself, creates a new GENRACEWAY_LIMITED or raises the existing
%      singleton*.
%
%      H = GENRACEWAY_LIMITED returns the handle to a new GENRACEWAY_LIMITED or the handle to
%      the existing singleton*.
%
%      GENRACEWAY_LIMITED('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in GENRACEWAY_LIMITED.M with the given input arguments.
%
%      GENRACEWAY_LIMITED('Property','Value',...) creates a new GENRACEWAY_LIMITED or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before GenRaceway_limited_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to GenRaceway_limited_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help GenRaceway_limited

% Last Modified by GUIDE v2.5 23-Apr-2019 22:39:02

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
    'gui_Singleton',  gui_Singleton, ...
    'gui_OpeningFcn', @GenRaceway_limited_OpeningFcn, ...
    'gui_OutputFcn',  @GenRaceway_limited_OutputFcn, ...
    'gui_LayoutFcn',  [] , ...
    'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
function GenRaceway_limited_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to GenRaceway_limited (see VARARGIN)

% Choose default command line output for GenRaceway_limited
handles.output = hObject;

%update Ax and introduce mousedown callback to it. Initalize Fig userdata
%for plotedit toolbar
Ax=handles.Ax;
grid(Ax,'on'); hold(Ax,'on'); axis(Ax,'manual');

%define GuiStruct and input it into handles
GuiStruct.CarCGH=gobjects(1); delete(GuiStruct.CarCGH);
GuiStruct.CarArrowH=gobjects(1); delete(GuiStruct.CarCGH);
GuiStruct.CarCG=[]; %[X,Y] of car
GuiStruct.CarN=[]; %[X,Y] normalized direction of car
GuiStruct.UnOrderedYellowCones=[]; %[X,Y] matrix
GuiStruct.UnOrderedBlueCones=[]; %[X,Y] matrix
GuiStruct.OrderedYellowCones=[]; %[X,Y] matrix
GuiStruct.OrderedBlueCones=[]; %[X,Y] matrix

handles.GuiStruct=GuiStruct; %input into handles
guidata(hObject, handles); %save handles in figure
function varargout = GenRaceway_limited_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;
%% Callbacks
function RandPush_Callback(hObject, eventdata, handles)
handles=guidata(hObject); %Obtain updated handles
Ax=handles.Ax; Yellow=[1,1,0]; Blue=[0,0,1];

%Yellow cones:
YSct=findobj(Ax,'userdata',Yellow); %find yellow cones by userdata
Yamnt=length(YSct); %obtain amount of yellow cones
Yxy=cell2mat(arrayfun(@(x) [x.XData,x.YData],YSct,'un',0)); %obtain XY of yellow cones
Yorder=randperm(Yamnt); %create random indices vector
YSct=YSct(Yorder); %reorder yellow cones
Ycolors=autumn(Yamnt); %obtain colors of cones by colormap
for k=1:Yamnt %recolor cones to fit new order
    YSct(k).CData=Ycolors(k,:);
end

%Blue cones:
BSct=findobj(Ax,'userdata',Blue); %find blue cones by userdata
Bamnt=length(BSct); %obtain amount of blue cones
Bxy=cell2mat(arrayfun(@(x) [x.XData,x.YData],BSct,'un',0)); %obtain XY of blue cones
Border=randperm(Bamnt); %create random indices vector
BSct=BSct(Border); %reorder blue cones
Bcolors=winter(Bamnt); %obtain colors of cones by colormap
for k=1:Bamnt %recolor cones to fit new order
    BSct(k).CData=Bcolors(k,:);
end

%Input data into GuiStruct
handles.GuiStruct.UnOrderedYellowCones=Yxy;
handles.GuiStruct.UnOrderedBlueCones=Bxy;

%allow for Order cones push to activate
Teal=[0.7,0.8,1];
handles.OrderConesPush.BackgroundColor=Teal;

guidata(hObject, handles); %save handles in figure
function OrderConesPush_Callback(hObject, eventdata, handles)
handles=guidata(hObject); %Obtain updated handles
Red=[1,0.6,0.6];
if norm(hObject.BackgroundColor-Red)<eps
    errordlg('Please Randomize the cones first','A-lon'); return
end
if ~isvalid(handles.GuiStruct.CarCGH) %check if car is valid
    errordlg('Please place a car first','A-lon'); return
end

Ax=handles.Ax;
CarCG=handles.GuiStruct.CarCG; %[X,Y]
CarN=handles.GuiStruct.CarN; %[X,Y] (direction)
Yxy=handles.GuiStruct.UnOrderedYellowCones;
Bxy=handles.GuiStruct.UnOrderedBlueCones;

%method returns indcies of relevant cones by order
switch handles.MethodPopMenu.Value
    case 1 %Traveling Salesman
        errordlg('Method Undeveloped','A-lon');
    case 2 %Delaunay Triangulation
        [Yind,Bind]=DelaunayTriangulation(Ax,Yxy,Bxy,CarCG,CarN);
    case 3 %Cone hopping
        errordlg('Method Undeveloped','A-lon');
        return
    case 4 %Double Arc fitting
        errordlg('Method Undeveloped','A-lon');
        return
end

%re-paint the cones to order
RePaintCones(Ax,Yind,Bind)

%save ordered cones into handles to be interpolated later
handles.GuiStruct.OrderedYellowCones=Yxy(Yind,:);
handles.GuiStruct.OrderedBlueCones=Bxy(Bind,:);

%allow for Intrp push to activate
Teal=[0.7,0.8,1];
handles.IntrpPush.BackgroundColor=Teal;
handles.ExportPush.BackgroundColor=Teal;

guidata(hObject, handles); %save handles in figure
function IntrpPush_Callback(hObject, eventdata, handles)
handles=guidata(hObject); %Obtain updated handles
Ax=handles.Ax; DarkYellow=0.9*[1,1,0]; Blue=[0,0,1];

%Obtain ordered cones xy
Yxy=handles.GuiStruct.OrderedYellowCones;
Bxy=handles.GuiStruct.OrderedBlueCones;

%Plot yellow
before = findall(Ax);
% if ~isempty(Yxy), fnplt(cscvn(Yxy')); end
if ~isempty(Yxy), plot(Yxy(:,1),Yxy(:,2)); end
after = setdiff(findall(Ax), before);
set(after,'color',DarkYellow,'lines','-','linew',2);

%Plot Blue
before = findall(Ax);
% if ~isempty(Bxy), fnplt(cscvn(Bxy')); end
if ~isempty(Bxy), plot(Bxy(:,1),Bxy(:,2)); end
after = setdiff(findall(Ax), before);
set(after,'color',Blue,'lines','-','linew',2);

%allow for Export push to activate
Teal=[0.7,0.8,1];
handles.ExportPush_Callback.BackgroundColor=Teal;
function ExportPush_Callback(hObject, eventdata, handles)
Red=[1,0.6,0.6];
if norm(hObject.BackgroundColor-Red)<eps
    errordlg('Please Randomize the cones first','A-lon'); return
end

handles=guidata(hObject); %Obtain updated handles
OYxy=handles.GuiStruct.OrderedYellowCones; assignin('base','OYxy',OYxy);
OBxy=handles.GuiStruct.OrderedBlueCones; assignin('base','OBxy',OBxy);
UnOYxy=handles.GuiStruct.UnOrderedYellowCones; assignin('base','UnOYxy',UnOYxy);
UnOBxy=handles.GuiStruct.UnOrderedBlueCones; assignin('base','UnOBxy',UnOBxy);
CarCG=handles.GuiStruct.CarCG; assignin('base','CarCG',CarCG);
CarN=handles.GuiStruct.CarN; assignin('base','CarN',CarN);
%mouse down callback
function Mdown(Ax,event,handles)
handles=guidata(Ax); %Obtain updated handles

if handles.PlaceConesToggle.Value %Place cones
    Color=handles.ColorPush.BackgroundColor;
    P=Ax.CurrentPoint; %returns 2x3 matrix of [x,y,z;x,y,z] in [front;back] format
    scatter(Ax,P(1,1),P(1,2),40,Color,'filled','linew',1,...
        'markeredgecolor','k','UserData',Color);
end
if handles.PlaceCarToggle.Value %Place car
    delete(handles.GuiStruct.CarCGH); %delete old car if existed
    delete(handles.GuiStruct.CarArrowH); %delete old quiver direction
    P=Ax.CurrentPoint; %returns 2x3 matrix of [x,y,z;x,y,z] in [front;back] format
    CarCG=[P(1,1),P(1,2)];
    CarN=[0,1];
    Purple=[0.7,0,0.7];
    L=0.15*mean([Ax.XLim(2)-Ax.XLim(1),Ax.YLim(2)-Ax.YLim(1)]);
    CarArrowH=quiver(Ax,CarCG(1),CarCG(2),L*CarN(1),L*CarN(2),...
        'linew',2,'color',Purple+0.1,'maxheadsize',0.6);
    CarCGH=scatter(Ax,P(1,1),P(1,2),50,Purple,'filled','linew',1,...
        'markeredgecolor','k');
    handles.GuiStruct.CarN=CarN; %reset car direction
    handles.GuiStruct.CarCG=CarCG;
    handles.GuiStruct.CarCGH=CarCGH;
    handles.GuiStruct.CarArrowH=CarArrowH;
end
if handles.RotateCarToggle.Value %Rotate Car
    P=Ax.CurrentPoint;
    CarCG=handles.GuiStruct.CarCG;
    CarN=RotateCar(CarCG,[P(1,1),P(1,2)]);
    L=0.15*mean([Ax.XLim(2)-Ax.XLim(1),Ax.YLim(2)-Ax.YLim(1)]);
    CarArrowH=handles.GuiStruct.CarArrowH;
    CarArrowH.UData=L*CarN(1); CarArrowH.VData=L*CarN(2);
    handles.GuiStruct.CarN=CarN;
end

guidata(Ax,handles); %save handles in figure
%simple callbacks
function PlaceCarToggle_Callback(hObject, eventdata, handles)
if hObject.Value == 1
    set(handles.Ax,'ButtonDownFcn',{@Mdown,handles})
    handles.Fig.Pointer='crosshair';
    handles.PlaceConesToggle.Value=0;
    handles.RotateCarToggle.Value=0;
    %disallow for some push buttons to activate
    Red=[1,0.6,0.6];
    handles.OrderConesPush.BackgroundColor=Red;
    handles.ExportPush.BackgroundColor=Red;
    handles.IntrpPush.BackgroundColor=Red;
else
    set(handles.Ax,'ButtonDownFcn','')
    handles.Fig.Pointer='arrow';
end
function RotateCarToggle_Callback(hObject, eventdata, handles)
handles=guidata(hObject); %Obtain updated handles
if ~isvalid(handles.GuiStruct.CarCGH) %if no car is placed
    errordlg('No car to rotate','A-lon');
    return
end
if hObject.Value == 1
    set(handles.Ax,'ButtonDownFcn',{@Mdown,handles})
    handles.Fig.Pointer='circle';
    handles.PlaceConesToggle.Value=0;
    handles.PlaceCarToggle.Value=0;
    %disallow for some push buttons to activate
    Red=[1,0.6,0.6];
    handles.OrderConesPush.BackgroundColor=Red;
    handles.ExportPush.BackgroundColor=Red;
    handles.IntrpPush.BackgroundColor=Red;
else
    set(handles.Ax,'ButtonDownFcn','')
    handles.Fig.Pointer='arrow';
end
function PlaceConesToggle_Callback(hObject, eventdata, handles)
if hObject.Value == 1
    set(handles.Ax,'ButtonDownFcn',{@Mdown,handles})
    handles.Fig.Pointer='hand';
    handles.PlaceCarToggle.Value=0;
    handles.RotateCarToggle.Value=0;
    %disallow for some push buttons to activate
    Red=[1,0.6,0.6];
    handles.OrderConesPush.BackgroundColor=Red;
    handles.ExportPush.BackgroundColor=Red;
    handles.IntrpPush.BackgroundColor=Red;
else
    set(handles.Ax,'ButtonDownFcn','')
    handles.Fig.Pointer='arrow';
end
function ClearPush_Callback(hObject, eventdata, handles)
cla(handles.Ax);
%disallow for some push buttons to activate
Red=[1,0.6,0.6];
handles.OrderConesPush.BackgroundColor=Red;
handles.ExportPush.BackgroundColor=Red;
handles.IntrpPush.BackgroundColor=Red;
function ResetColorsPush_Callback(hObject, eventdata, handles)
Ax=handles.Ax; Yellow=[1,1,0]; Blue=[0,0,1];

%Yellow cones:
YSct=findobj(Ax,'userdata',Yellow); %find yellow cones by userdata
Yamnt=length(YSct); %obtain amount of yellow cones
for k=1:Yamnt %recolor cones to fit new order
    YSct(k).CData=Yellow;
end

%Blue cones:
BSct=findobj(Ax,'userdata',Blue); %find blue cones by userdata
Bamnt=length(BSct); %obtain amount of blue cones
for k=1:Bamnt %recolor cones to fit new order
    BSct(k).CData=Blue;
end

%disallow for some push buttons to activate
Red=[1,0.6,0.6];
handles.OrderConesPush.BackgroundColor=Red;
handles.ExportPush.BackgroundColor=Red;
handles.IntrpPush.BackgroundColor=Red;
%deadsimple callbacks
function LoadConesPush_Callback(hObject, eventdata, handles)
handles=guidata(hObject); %Obtain updated handles

%Initalize
Ax=handles.Ax;
Cones=[]; %initalize    
uiload; %load Cones variable - mx3 [x,y,color] array
if isempty(Cones), return, end

%Extract XY from Cones
Bxy=Cones(Cones(:,3)==98,1:2);
Yxy=Cones(Cones(:,3)==121,1:2);


%Draw cones and store in handles
h=waitbar(0,'Please Wait. . . '); %introduce waitbar
axis(Ax,'auto'); %turn Ax limits to auto for drawing

%Draw yellow cones
handles.GuiStruct.UnOrderedYellowCones=Yxy(:,1:2);
Color=[1,1,0]; %yellow
for k=1:size(Yxy,1)
    scatter(Ax,Yxy(k,1),Yxy(k,2),40,Color,'filled','linew',1,...
        'markeredgecolor','k','UserData',Color);
end
waitbar(0.5,h); %updat waitbar

%draw blue cones
handles.GuiStruct.UnOrderedBlueCones=Bxy(:,1:2);
Color=[0,0,1]; %blue
for k=1:size(Bxy,1)
    scatter(Ax,Bxy(k,1),Bxy(k,2),40,Color,'filled','linew',1,...
        'markeredgecolor','k','UserData',Color);
end

delete(h); %kill waitbar
AutoLimsCheckBox_Callback(handles.AutoLimsCheckBox,eventdata,handles) %return Ax lims to what user decided
guidata(hObject, handles); %save handles in figure
function AutoLimsCheckBox_Callback(hObject, eventdata, handles)
Ax=handles.Ax;

if hObject.Value
    axis(Ax,'auto');
else
    axis(Ax,'manual');
end
function AxLimEdit_Callback(hObject, eventdata, handles)
Ax=handles.Ax;
lims=str2num(hObject.String);
if isempty(lims), errordlg('Please input a correct limit input','A-lon'); return, end
if length(lims)~=4,  errordlg('Please input a correct limit input','A-lon'); return, end

Ax.XLim=lims(1:2); Ax.YLim=lims(3:4);
function PlotEditToggle_Callback(hObject, eventdata, handles)
plotedit(handles.Fig,'on');
function ColorPush_Callback(hObject, eventdata, handles)
Yellow=[1,1,0]; Blue=[0,0,1];
if norm(hObject.BackgroundColor-Yellow)<eps %button is yellow, switch it blue
    hObject.BackgroundColor=Blue;
else %button is blue, make it yellow
    hObject.BackgroundColor=Yellow;
end
%% Method Functions
function [Yind,Bind]=DelaunayTriangulation(Ax,Yxy,Bxy,CarCG,CarVel)
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

%Initalize for output
Yind=[];
Bind=[];

%ask for user input
Prompt={'DrawLines (Bool)','Max Iterations','Cost Threshold'...
    ,'Radii Ratio Threshold','Pause Time','SphereFilter R','CarLength'};
Default={'0','10','-0.2','10','0','100000000','5'};
Uinput=inputdlg(Prompt,'Delaunay Triangulation Parameters',1,Default);
if isempty(Uinput), return, end
if any(cellfun(@isempty,Uinput)), errordlg('Please fill all inputs','A-lon'), return, end

%running parameters
PlotBool=str2double(Uinput{1});
MaxItrAmnt=str2double(Uinput{2});
CostThreshold=str2double(Uinput{3}); %if Cost>CostThreshold, exit
RRatioThreshold=str2double(Uinput{4}); %if Rbinding/Rbounded>RRatioThreshold, exit
PauseTime=str2double(Uinput{5});
SphereFilterR=str2double(Uinput{6});
CarLength=str2double(Uinput{7});

%initalize for functionatlity - passed Uinput "test"
Yind=nan(MaxItrAmnt,1);
Bind=nan(MaxItrAmnt,1);

%Preprocessing
CarDir=CarVel/norm(CarVel); %normalize velocity to unit direction
Yamnt=size(Yxy,1); Bamnt=size(Bxy,1);
YColor=121*ones(Yamnt,1); BColor=98*ones(Bamnt,1);
Cones=[[Bxy,BColor];[Yxy,YColor]]; %Create Cones-[x,y,color]
Cones=SphereFilter(Cones,CarCG,SphereFilterR);
%create color array that corelates with DT.Points
ConesColors=Cones(:,3); %Optimaly would sit with DelaunayTriangulation points. This is a compremise.

%IMPORTANT: DT.Points - exactly the same as Cones(:,[1,2])

%Triangulate
DT=delaunayTriangulation(Cones(:,[1,2])); %Triangulate only /w Cones
ID=pointLocation(DT,CarCG); %attempt to find triangle which contains CarCG
if isnan(ID) %if CarCG isoutside of convex hull
    DT=delaunayTriangulation([Cones(:,[1,2]);CarCG]); %Triangualte with CarCG as node
    ID=pointLocation(DT,CarCG+0.5*CarLength*CarDir); %find first triangle in Car direction
    ConesColors=[ConesColors;0]; %add 0 color for CarCG node
end
[NewID,Newu,NewCrossEdge]=TriOne(DT,ConesColors,ID,CarDir);
if isempty(NewID), return, end %crossed into no-man's land

%insert first cones into lists (Bind/Yind)
if ConesColors(NewCrossEdge(1))=='y', Yind(1)=NewCrossEdge(1);
else, Bind(1)=NewCrossEdge(1); end
if ConesColors(NewCrossEdge(2))=='y', Yind(1)=NewCrossEdge(2);
else, Bind(1)=NewCrossEdge(2); end
NewV=setdiff(DT.ConnectivityList(NewID,:),NewCrossEdge);
if ConesColors(NewV)=='y', Yind(2)=NewV;
else, Bind(2)=NewV; end

%plot triangulation and first triangle if user decided
if PlotBool
    %all triangulation
    triplot(DT,'color','k','tag','triplot','deletefcn',{@deleteDT,Ax});
    %first triangle
    V=DT.ConnectivityList(ID,:);
    triplot(V,DT.Points(:,1),DT.Points(:,2),'color','r','linew',2,'tag','triplot',...
        'deletefcn',{@deleteDT,Ax});
    quiver(CarCG(1),CarCG(2),+0.15*CarDir(1),+0.15*CarDir(2),'color',[1,0.8,0.2],...
        'linew',2,'MaxHeadSize',0.5,'tag','quiver','deletefcn',{@deleteDT,Ax})
end

for Itr=3:MaxItrAmnt
    if PlotBool %if user wants to plot - plot it all
        V=DT.ConnectivityList(NewID,:);
        triplot(V,DT.Points(:,1),DT.Points(:,2),'color','r','linew',2,'tag','triplot',...
            'deletefcn',{@deleteDT,Ax});
        p=mean(DT.Points(NewCrossEdge,:));
        quiver(p(1)-0.075*Newu(1),p(2)-0.075*Newu(2),+0.15*Newu(1),+0.15*Newu(2),'color',[1,0.8,0.2],...
            'linew',2,'MaxHeadSize',0.5,'tag','quiver','deletefcn',{@deleteDT,Ax})
        pause(PauseTime);
    end
    
    %find next triangle
    [NewID,Newu,NewCrossEdge,Cost,RRatio]=FindNextTriangle(DT,ConesColors,NewID,Newu,NewCrossEdge);
    %check conditions, if not good enough - break
    if isempty(NewID) || CostThreshold<Cost || RRatioThreshold<RRatio, break; end
    
    %add next cone to Bind/Yind
    NewV=setdiff(DT.ConnectivityList(NewID,:),NewCrossEdge);
    if ConesColors(NewV)=='y', Yind(Itr)=NewV;
    else, Bind(Itr)=NewV; end
end

if PlotBool && Itr==MaxItrAmnt %if user wants to plot - plot last triangle
    V=DT.ConnectivityList(NewID,:);
    triplot(V,DT.Points(:,1),DT.Points(:,2),'color','r','linew',2,'tag','triplot',...
        'deletefcn',{@deleteDT,Ax});
    p=mean(DT.Points(NewCrossEdge,:));
    quiver(p(1)-0.075*Newu(1),p(2)-0.075*Newu(2),+0.15*Newu(1),+0.15*Newu(2),'color',[1,0.8,0.2],...
        'linew',2,'MaxHeadSize',0.5,'tag','quiver','deletefcn',{@deleteDT,Ax})
    pause(PauseTime);
end

Bind(any(isnan(Bind),2),:)=[]; Yind(any(isnan(Yind),2),:)=[]; %delete rows with NaNs
Bamnt=sum(Cones(:,3)==98); %for drawing - might have changed because of SphereFilter
Yind=Yind-Bamnt; %Bamnt was recalculated if hemisphere filter was activated
Yind=Yind'; Bind=Bind'; %make them row vectors - expected output for RepaintCones
%% Some more functions
function RePaintCones(Ax,Yind,Bind)
%Yind, Bind - row vectors of indcies pointing on cones by order
Yellow=[1,1,0]; Blue=[0,0,1]; Gray=0.5*[1,1,1];
YSct=findobj(Ax,'userdata',Yellow); %find yellow cones by userdata
BSct=findobj(Ax,'userdata',Blue); %find blue cones by userdata

Bamnt=length(Bind);
Bcolors=winter(Bamnt); %obtain colors of cones by colormap
i=1; %for color order
for k=Bind %recolor cones to fit new order
    BSct(k).CData=Bcolors(i,:);
    i=i+1;
end
Bdiff=setdiff(1:length(BSct),Bind); %set cones that arent in Bind to gray
if ~isempty(Bdiff)
    set(BSct(Bdiff),'CData',Gray);
end

Yamnt=length(Yind);
Ycolors=autumn(Yamnt); %obtain colors of cones by colormap
i=1;
for k=Yind %recolor cones to fit new order
    YSct(k).CData=Ycolors(i,:);
    i=i+1;
end
Ydiff=setdiff(1:length(YSct),Yind); %set cones that arent in Bind to gray
if ~isempty(Ydiff)
    set(YSct(Ydiff),'CData',Gray);
end
function CarN=RotateCar(CarCG,P)
%Input:
%CarCg - [X,Y]
%P - [X,Y] where user clicked

% https://de.mathworks.com/matlabcentral/answers/180131-how-can-i-find-the-angle-between-two-vectors-including-directional-information
% old formulation for transformation matrix - Calcualte CarN
% v=[0,1];
% Theta=atan2(-CarN(1)*v(2)+CarN(2)*v(1),CarN(1)*v(1)+CarN(2)*v(2)); %atan(sin(theta)/cos(theta).
% R=makehgtform('zRotate',Theta);
% T=makehgtform('Translate',[CarCG(1),CarCG(2),0]);
% CarTransfm.Matrix=T*R;
CarN=(P-CarCG)/norm(P-CarCG);
%for Traveling Salesman
function [ConesXY,I]=SortConesbyCar(ConesXY,CarXY)
%Sorts an already orginized array of cones by closest cone to car
%returns indcies of sorted ConesXY
PrevD=norm(ConesXY(1,:)-CarXY);
Ind=1;
for k=2:size(ConesXY,1)
    D=norm(ConesXY(k,:)-CarXY);
    if D<PrevD
        PrevD=D;
        Ind=k;
    end
end
ConesXY=circshift(ConesXY,-Ind,1);
I=circshift(1:size(ConesXY,1),-Ind,2);
%for Delaunay Triangulation
function [NewID,Newu,NewCrossEdge]=TriOne(DT,ConesColors,ID,CarDir)
%check the cost of three 

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
%NewV - %New Vertex in new triangle (=new triangle V/NewCrossEdge V).

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
function J=EdgeCost(InCenter,P,u)
%InCenter - in center of triangle [x,y]
%P=[x,y,color] of edges vertices 2x3
%u - direction from previous triangle to current triangle
J=0;
if abs((P(1,3)-P(2,3)))<eps %if they are the same color
    J=J+1; else, J=J-1; end
J=J-dot(u,OutFacingNormal(P(:,1:2),InCenter)); %+bad points for difference in direction
function n=OutFacingNormal(P,InCenter)
%function returns the out facing normal for an edge of a triangle [p1,p2]
%with InCenter.
%P=[x,y] of edge's vertices 2x2
q=PorjPnt2Line(P(1,:),P(2,:),InCenter);
n=(q-InCenter)/norm(q-InCenter);
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
function projq=PorjPnt2Line(p1,p2,q)
%p1,p2,q,projq - [X,Y] points.
%p1,p2 - represent a line
t=(p2-p1)/norm(p2-p1);
projq=dot((q-p1),t)*t+p1;
function deleteDT(hObject,~,Ax)
%if one quiver/triplot is deleted, delete all of them
delete([findobj(Ax,'Tag','quiver');findobj(Ax,'Tag','triplot')]);
%% empty GUI functions
function AxLimEdit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to AxLimEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function RandPauseEdit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to RandPauseEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function MethodPopMenu_Callback(hObject, eventdata, handles)
function MethodPopMenu_CreateFcn(hObject, eventdata, handles)
% hObject    handle to MethodPopMenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
%% functions not in use
function [FilteredCones,FInd]=HemiSphereFilter(Cones,CarCG,CarDir,R)
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
FInd=(SqDistance<R^2) & ConesR*CarDir'>0;
FilteredCones=Cones(FInd,:);
function [NewID,Newu,NewCrossEdge]=TriOne_CarOut(DT,ConesColors,ID)
%Major assumption:
%Car is built of yellow cone - blue cone- car
%finds yellow and blue cone and declares as edge to cross

%Input:
%DT - DelaunayTriangulation containing DT.Points and DT.ConnectivityList
%ConesColors - colors of cones (98 for blue, 121 for yellow) with same
%indexing as DT.Points
%ID - number of triangle in DT.ConnectivityList (row index)

%Output:
%NewID - number of new triangle in DT.ConnectivityList (row index)
%Newu - direction of enterance to new triangle (NewID)
%NewCrossEdge - [V1,V2] of edge that we crossed to get from ID->NewID
%V1 and V2 refer to vertex indcies (row) in DT.Points

V=DT.ConnectivityList(ID,:); %find vertex indcies of ID
ConesColors=[ConesColors;0]; %add car node color as 0
BCone=V(ConesColors(V)==98); YCone=V(ConesColors(V)==121); %find blue and yellow cone in triangle
NewCrossEdge=[BCone,YCone]; %
InCenter=incenter(DT,ID); %find incenter of ID
Newu=OutFacingNormal(DT.Points(NewCrossEdge,:),InCenter);
Attchments=edgeAttachments(DT,NewCrossEdge); %cell containing 2 triangles ID {current,next}
NewID=setdiff(Attchments{1},ID);
