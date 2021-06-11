function OMxy=FindMiddle(OBxy,OYxy)
%OBxy,OYxy,OMxy - same format
%OBxy - ordered blue cones.
%OYxy - ordered yellow cones.
%OMxy - ordered middle points
%[x,y] size mx2 nuermic matrices containing coordiantes of cones

%decide on inner and outer cones by amount (infantile)
Bamnt=size(OBxy,1); Yamnt=size(OYxy,1);
[Mamnt,Ind]=min([Bamnt,Yamnt]); %amount of middle points - as inner track cones
if Ind==1, Inxy=OBxy; Outxy=OYxy;
else, Inxy=OYxy; Outxy=OBxy; end

OMxy=zeros(Mamnt,2);
for k=1:Mamnt
    InCone=Inxy(k,:);
    [~,Ind]=min(diag((Outxy-InCone)*(Outxy-InCone)')); %find closest OutCone to k-th InCone
    OutCone=Outxy(Ind,:);
    OMxy(k,:)=0.5*(InCone+OutCone);
end
end