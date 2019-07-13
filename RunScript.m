%% Peaks
PtchAmnt=16;
Layers=sqrt(PtchAmnt);
BezO=9; 
N=Layers*BezO+1;
[X,Y,Z]=peaks(N);
MeshNodes=zeros(size(X));
MeshNodes(:,:,1)=X; MeshNodes(:,:,2)=Y; MeshNodes(:,:,3)=Z;
CP=BezCP(MeshNodes,BezO,'Method','Block');
PseudoInverseCP=CP.PesudoInverseVertices;

%compare PeaksCP to original surface
%PeaksCP - mesh nodes of peaks surface are the control points for the bezier surface mesh
Ax11=CP.DrawBezierPatches('color',[1,0,0],'facealpha',1,'title','CP Bezier mesh vs Peaks Surface');
BezCP.DrawPointCloud(CP.Vertices,'color',[0,1,0],'msize',20,'Ax',Ax11); %draw control points
surf(Ax11,X,Y,Z,'EdgeColor','none','FaceAlpha',0.5)

%compare PseudoPeaksCP to original peaks surface
%PseudoPeaksCP - bezier surface mesh attempts to equal peaks surface
Ax12=PseudoInverseCP.DrawBezierPatches('color',[1,1,1],'title','PsuedoInverseCP Bezier mesh vs Peaks Surface');
BezCP.DrawPointCloud(PseudoInverseCP.Vertices,'color',[0,1,0],'msize',20,'Ax',Ax12); %draw control points
surf(Ax12,X,Y,Z,'EdgeColor','none','FaceAlpha',1)

Ax13=PseudoInverseCP.DrawBezierPatches('color',[1,1,1],'title','PsuedoInverseCP Bezier mesh vs CP Bezier mesh');
CP.DrawBezierPatches('color',[1,0,0],'facealpha',1,'Ax',Ax13);

%obtain figure handles of recently created axes
h11=Ax11.Parent; h12=Ax12.Parent; h13=Ax13.Parent;

%move axes into subplot array
fig=figure('color',[0,0,0],'units','normalized','outerposition',[0 0 1 1]);
subplot(1,3,1,Ax11,'parent',fig);
subplot(1,3,2,Ax12,'parent',fig);
subplot(1,3,3,Ax13,'parent',fig);

%delete old figures
close(h11,h12,h13);
clear('Ax11','Ax12','Ax13');

%% Calculate 
Stmp=Bez4Stmp('Roie.stl','Cap',true,'SphLayers',1,'CylLayers',1,'Slices',4,'BezierOrder',3,...
    'XcenterCalculationMethod','normalSTD');

%% Update
Stmp.Cap=0;
Stmp.Slices=10;
Stmp.SphLayers=3;
Stmp.CylLayers=2;
Stmp.BezierOrder=3;
Stmp=Stmp.UpdateObj;

%% Draw Patches
pausetime=0;
% Stmp.CP.DrawBezierPatches('pausetime',pausetime);
% Stmp.PsuedoInverseCP.DrawBezierPatches('pausetime',pausetime);
Stmp.PsuedoInverseCP.DrawBezierPatches('curvature','gaussian','PauseTime',0);

%% Draw point clouds - orginial, compact and surfaces
Ax=Stmp.DrawPointCloud(Stmp.PointCloud,'color',[1,0,0],'msize',20); %original
Stmp.DrawPointCloud(Stmp.Compact,'color',[0,0,1],'msize',20,'Ax',Ax); %compact
Stmp.DrawPointCloud(Stmp.PsuedoInverseCP.CombinePatches(30),'color',[1,1,1],'Ax',Ax); %compact
Stmp.DrawPointCloud(Stmp.PsuedoInverseCP.Vertices,'color',[0,1,0],'Ax',Ax,'msize',20); %CP vertices

%% Hausdorff distance with plot
P=Stmp.PointCloud.Location;
Q=Stmp.PsuedoInverseCP.CombinePatches(30);
szP=size(P); szQ=size(Q);
if numel(szP)==3, P=reshape(P,szP(1)*szP(2),3); end
if numel(szQ)==3, Q=reshape(Q,szQ(1)*szQ(2),3); end
Threshold=30;
P=P(P(:,3)>Threshold,:); Q=Q(Q(:,3)>Threshold,:); %filter buttom noise
[hd,pInd,qInd]=Stmp.Hausdorff(P,Q); %<----------hausdorff distance
Phd=P(pInd,:); Qhd=Q(qInd,:);
%calculate radial distance (radial off mean point of Phd and Qhd)
m=mean([Phd;Qhd]);
Xcntr=Stmp.Xcenter;
if m(3)>Xcntr(3), t=(m-Xcntr)/norm(m-Xcntr,2);
else, t=[(m(1:2)-Xcntr(1:2)),0]/norm(m(1:2)-Xcntr(1:2),2); end
rhd=abs(dot(Phd-Qhd,t));
%plot
Ax=Stmp.DrawPointCloud(Stmp.PointCloud,'color',[0,0,1],'msize',15); %original
Stmp.DrawPointCloud(Stmp.PsuedoInverseCP.CombinePatches(30),'color',[1,1,1],'Ax',Ax); %compact
Stmp.DrawPointCloud([Phd;Qhd],'color',[1,0,0],'msize',20,'Ax',Ax,...
    'title',sprintf('Hausdorff distance %.2g with Radial displacement of %.2g',hd,rhd)); %compact