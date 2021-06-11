%% Testing OrderConesDT

%GUI steps:
%Step1: Place cones and car in GenRaceway
%Step2: Solve with Delanuay Triangulation
%Step3: Use Export2WS to Obtain variables
%OBxy,OYxy,UnOBxy,UnOYxy,CarCG,CarN

%preprocessing
Yamnt=size(UnOYxy,1); Bamnt=size(UnOBxy,1);
YColor=121*ones(Yamnt,1); BColor=98*ones(Bamnt,1);
Cones=[[UnOBxy,BColor];[UnOYxy,YColor]]; %Create Cones-[x,y,color]

%Apply OrderConesDT
[BCones,YCones]=OrderConesDT(Cones,CarCG,CarN);

%Plot
Fig=figure;
AxTrack=axes(Fig);
hold(AxTrack,'on'); grid(AxTrack,'on'); axis(AxTrack,'manual');
scatter(AxTrack,BCones(:,1),BCones(:,2),40,[0,0,1],'filled','linew',1,...
    'markeredgecolor','k');
scatter(AxTrack,YCones(:,1),YCones(:,2),40,[1,1,0],'filled','linew',1,...
    'markeredgecolor','k');
plot(AxTrack,YCones(:,1),YCones(:,2),'color',[1,1,0],'linew',2);
plot(AxTrack,BCones(:,1),BCones(:,2),'color',[0,0,1],'linew',2);

%Compare with GUI image
%% Testing BuildTrack

%testing data
OBxy=[0.3386    0.3149 %13 cones total
    0.4452    0.2046
    0.6009    0.1619
    0.7507    0.2117
    0.8112    0.3470
    0.8112    0.4715
    0.7622    0.6246
    0.6671    0.6957
    0.5317    0.7064
    0.3876    0.6851
    0.2983    0.5676
    0.2550    0.4253
    0.3386    0.3149];

OYxy=[  0.2233    0.1263 %20 cones total
    0.2896    0.0552
    0.4251    0.0374
    0.5461    0.0374
    0.7190    0.0765
    0.8919    0.1192
    0.9352    0.2260
    0.9496    0.3577
    0.9409    0.5534
    0.9352    0.7313
    0.8545    0.8843
    0.7190    0.9306
    0.5403    0.9448
    0.3242    0.9235
    0.2205    0.8665
    0.1225    0.7384
    0.0908    0.5996
    0.1052    0.4466
    0.1455    0.3185
    0.1628    0.2544];

%Plot ground truth
FigGT=figure;
AxGT=axes(FigGT); %ground truth
hold(AxGT,'on'); grid(AxGT,'on'); axis(AxGT,'manual');
scatter(AxGT,OBxy(:,1),OBxy(:,2),40,[0,0,1],'filled','linew',1,...
    'markeredgecolor','k');
scatter(AxGT,OYxy(:,1),OYxy(:,2),40,[1,1,0],'filled','linew',1,...
    'markeredgecolor','k');

%split to groups. %13 blue cones total, and 20 yellow cones total
BFrames{1}=OBxy(1:4,:);
BFrames{2}=OBxy(3:7,:);
BFrames{3}=OBxy(5:10,:);
BFrames{4}=OBxy(10:13,:);

YFrames{1}=OYxy(1:7,:);
YFrames{2}=OYxy(5:12,:);
YFrames{3}=OYxy(8:14,:);
YFrames{4}=OYxy(14:18,:);
YFrames{5}=OYxy(18:20,:);

%implament BuildTrack
% BTrackCones=BuildTrack(BFrames);
% YTrackCones=BuildTrack(YFrames);

%plot built track
FigT=figure;
AxT=axes(FigT); %ground truth
hold(AxT,'on'); grid(AxT,'on'); axis(AxT,'manual');
scatter(AxT,BTrackCones(:,1),BTrackCones(:,2),40,[0,0,1],'filled','linew',1,...
    'markeredgecolor','k');
scatter(AxT,YTrackCones(:,1),YTrackCones(:,2),40,[1,1,0],'filled','linew',1,...
    'markeredgecolor','k');
%% Trying algorithms from simulator and MidPoints

load('SimAllCones.mat');
CarCG=1e4*[-1.7540,-1.3050];
CarN=[0.9944,-0.1052];

%Apply OrderConesDT
[BCones,YCones]=OrderConesDT(Cones,CarCG,CarN);
MidPoints=FindMiddle(BCones,YCones);

%Plot cones
AxCones=subplot(2,1,1);
hold(AxCones,'on'); grid(AxCones,'on');
scatter(AxCones,Bxy(:,1),Bxy(:,2),40,[0,0,1],'filled','linew',1,...
    'markeredgecolor','k');
scatter(AxCones,Yxy(:,1),Yxy(:,2),40,[1,1,0],'filled','linew',1,...
    'markeredgecolor','k');
scatter(AxCones,CarCG(1),CarCG(2),40,[1,0,1],'filled','linew',1,...
    'markeredgecolor','k');
quiver(CarCG(1),CarCG(2),2000*CarN(1),2000*CarN(2),'linew',2,'color',[0.8,0.8,0]);
%Plot track
AxTrack=subplot(2,1,2);
hold(AxTrack,'on'); grid(AxTrack,'on');
scatter(AxTrack,CarCG(1),CarCG(2),40,[1,0,1],'filled','linew',1,...
    'markeredgecolor','k');
quiver(CarCG(1),CarCG(2),2000*CarN(1),2000*CarN(2),'linew',2,'color',[0.8,0.8,0]);
plot(AxTrack,YCones(:,1),YCones(:,2),'color',[1,1,0],'linew',2);
plot(AxTrack,BCones(:,1),BCones(:,2),'color',[0,0,1],'linew',2);
plot(AxTrack,MidPoints(:,1),MidPoints(:,2),'color',[0.5,0,0],'linew',2);
%% testing with SimFrames

load('SimFrames.mat');
AxFrames=subplot(2,1,1); hold(AxFrames,'on'); grid(AxFrames,'on');
for k=1:length(Frames)
[BFrames{k},YFrames{k}]=OrderConesDT(Frames{k},CarCG{k},CarN{k});
plot(AxFrames,YFrames{k}(:,1),YFrames{k}(:,2),'color',[1,1,0],'linew',0.5);
plot(AxFrames,BFrames{k}(:,1),BFrames{k}(:,2),'color',[0,0,1],'linew',0.5);
scatter(AxFrames,Frames{k}(:,1),Frames{k}(:,2),20,'filled')
end

%implament BuildTrack
BTrackCones=BuildTrack(BFrames(1:end));
YTrackCones=BuildTrack(YFrames(1:end));
MidPoints=FindMiddle(BTrackCones,YTrackCones); %find middle points

%Plot track
AxTrack=subplot(2,1,2);
hold(AxTrack,'on'); grid(AxTrack,'on');
plot(AxTrack,YTrackCones(:,1),YTrackCones(:,2),'color',[1,1,0],'linew',2);
plot(AxTrack,BTrackCones(:,1),BTrackCones(:,2),'color',[0,0,1],'linew',2);
plot(AxTrack,MidPoints(:,1),MidPoints(:,2),'color',[0.5,0,0],'linew',2);