%https://www.mathworks.com/help/robotics/ref/robotics.node.html

%Architecure
%node1(OrbSlam+Yolo) ---> node2(Order Cones)


%every tick:
%node1 sends unordered cones to node2 
%node2 orders the cones and saves in workspace. plotting also occurs for
%control

global Bxy Yxy CG N

%build architecture
master=robotics.ros.Core; %create ROS master

node1=robotics.ros.Node('/OrbSlam'); %create node 1
node2=robotics.ros.Node('/OrderCones'); %create node 2

pub=robotics.ros.Publisher(node1,'/UOcones','std_msgs/Float64MultiArray');
sub=robotics.ros.Subscriber(node2,'/UOcones','std_msgs/Float64MultiArray',@ROSSubCB);
%% implement sequence and plot it
load('Data/SimFrames');
Ax=axes;
axis(Ax,'manual'); xlim(Ax,[-2.2,-1.4]*10^4); ylim(Ax,[-3,-1.2]*10^4);
hold(Ax,'on'); grid(Ax,'on');
[BFrames,YFrames]=deal(cell(length(Frames),1));

%Initailize video file
MovieName='ShowCase';
vidWriter = VideoWriter([MovieName,'.mp4'],'MPEG-4');
vidWriter.FrameRate=2;
open(vidWriter);

for k=1:length(Frames)
    %Compress data
    ArrayFrame=Frames{k}(:);
    DataArray=[CarCG{k}(:);CarN{k}(:);ArrayFrame];
    %Send msg
    Msg=rosmessage(pub);
    Msg.Data=DataArray;
    send(pub,Msg)% ROSSubCB callback activated
    pause(0.1); %give it some time
    
    BFrames{k}=Bxy; YFrames{k}=Yxy;
    %plot
    plot(Ax,Bxy(:,1),Bxy(:,2),'color',[0,0,1],'linew',2);
    plot(Ax,Yxy(:,1),Yxy(:,2),'color',[1,1,0],'linew',2);
    scatter(Ax,Bxy(:,1),Bxy(:,2),40,[0,0,1],'filled','linew',1,...
        'markeredgecolor','k');
    scatter(Ax,Yxy(:,1),Yxy(:,2),40,[1,1,0],'filled','linew',1,...
        'markeredgecolor','k');
    scatter(Ax,CG(1),CG(2),40,[1,0,1],'filled','linew',1,...
        'markeredgecolor','k');
    quiver(CG(1),CG(2),2000*N(1),2000*N(2),'linew',2,'color',[0.8,0.8,0]);
    
    pause(0.1);
    writeVideo(vidWriter,getframe(Ax));
end

BTrackCones=BuildTrack(BFrames);
YTrackCones=BuildTrack(YFrames);
MidPoints=FindMiddle(BTrackCones,YTrackCones); %find middle point
plot(Ax,MidPoints(:,1),MidPoints(:,2),'color',[0.5,0,0],'linew',2);

writeVideo(vidWriter,getframe(Ax));
close(vidWriter);