function ROSSubCB(sub,msg)
global Bxy Yxy CG N
%Extract data
DataArray=msg.Data;
CG=DataArray(1:2)'; N=DataArray(3:4)';
Cones=reshape(DataArray(5:end),[(length(DataArray)-4)/3,3]);

%activate algorithm
[Bxy,Yxy]=OrderConesDT(Cones,CG,N);
end