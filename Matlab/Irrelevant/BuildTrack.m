function TrackCones=BuildTrack(Frames)
%IMPORTANT NOTE:
%this function does not distinguish color.
%to sort yellow and blue cones, one must call this function twice

%Input:
%Frames - cell array m X 1 containing km X 2 double arrays of cones
%positions

%Output:
%TrackCones - k X 2 double arrays of cones

R=eps; %cones width [m]
% TrackCones=nan(size(cell2mat(Frames(:)),1),1); %preallocate

Frames=Frames(:); %ensure column
FramesAmnt=length(Frames);
TrackCones=Frames{1}; %first frame inserted into output
for m=2:FramesAmnt
    FrmCone1=Frames{m}(1,:);
    Dist=diag((TrackCones-FrmCone1)*(TrackCones-FrmCone1)');
    Ind=find(Dist<R^2);
    if ~isempty(Ind), TrackCones=TrackCones(1:(Ind-1),:); end
    TrackCones=[TrackCones;Frames{m}];
end
end