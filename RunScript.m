%% Peaks
PtchAmnt=16;
Layers=sqrt(PtchAmnt);
BezO=4; 
N=Layers*BezO+1;
[X,Y,Z]=peaks(N);
MeshNodes=zeros(size(X));
MeshNodes(:,:,1)=X; MeshNodes(:,:,2)=Y; MeshNodes(:,:,3)=Z;
CP=BezCP(MeshNodes,BezO,'Method','Block');
PseudoInverseCP=CP.PesudoInverseVertices;

%move axes into subplot array
fig=figure('color',[0,0,0],'units','normalized','outerposition',[0 0 1 1]);
Ax11=subplot(1,2,1,BezCP.CreateDrawingAxes(fig),'parent',fig);
Ax12=subplot(1,2,2,BezCP.CreateDrawingAxes(fig),'parent',fig);

%compare PeaksCP to original surface
%PeaksCP - mesh nodes of peaks surface are the control points for the bezier surface mesh
srfH=CP.DrawBezierPatches('Ax',Ax11,'color',[0,0.7,0.7],'facealpha',1,...
    'edgecolor',0.5*[1,1,1],'title','Control Points of CP are vertices in Peaks Surface');
pcH=BezCP.DrawPointCloud(MeshNodes,'Ax',Ax11,'color',[1,0,0],'msize',20);
subset=[srfH(1),pcH];
legend(subset,'\color{white}CP Bezier mesh','\color{white}Peaks surface vertices & CP control points')

%compare PseudoPeaksCP to original peaks surface
%PseudoPeaksCP - bezier surface mesh attempts to equal peaks surface
srfH=PseudoInverseCP.DrawBezierPatches('Ax',Ax12,'color',[1,1,1],'facealpha',1,...
    'edgecolor',0.5*[1,1,1],'title','Control points of PsuedoInverseCP are outside Peaks Surface');
pcH1=BezCP.DrawPointCloud(PseudoInverseCP.Vertices,'Ax',Ax12,'color',[0,1,0],'msize',20); %draw control points
pcH2=BezCP.DrawPointCloud(MeshNodes,'Ax',Ax12,'color',[1,0,0],'msize',20);
subset=[srfH(1),pcH1,pcH2];
legend(subset,'\color{white}PseudoInverseCP Bezier mesh','\color{white}PsuedoInverseCP control points',...
    '\color{white}Peaks surface vertices')
%% Calculate 
Stmp=Bez4Stmp('Roie.stl','Cap',true,'SphLayers',2,'CylLayers',2,'Slices',4,'BezierOrder',3,...
    'XcenterCalculationMethod','normalSTD');

%% Update
Stmp.Cap=0;
Stmp.Slices=4;
Stmp.SphLayers=2;
Stmp.CylLayers=2;
Stmp.BezierOrder=3;
Stmp=Stmp.UpdateObj;

%% Draw Patches
fig=figure('color',[0,0,0]);
Ax=Stmp.CreateDrawingAxes(fig);

pausetime=0;
% Stmp.CP.DrawBezierPatches('pausetime',pausetime);
Stmp.PsuedoInverseCP.DrawBezierPatches('Ax',Ax,'pausetime',pausetime);
% Stmp.PsuedoInverseCP.DrawControlPoints('Ax',Ax,'pausetime',pausetime);
% Stmp.PsuedoInverseCP.DrawBezierPatches('curvature','gaussian','PauseTime',0);

%% Draw point clouds - orginial, compact and surfaces
fig=figure('color',[0,0,0]);
Ax=Stmp.CreateDrawingAxes(fig);
Stmp.DrawPointCloud(Stmp.PointCloud,'color',[1,0,0],'msize',20,'Ax',Ax); %original
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
fig=figure('color',[0,0,0]);
Ax=Stmp.CreateDrawingAxes(fig);
Stmp.DrawPointCloud(Stmp.PointCloud,'color',[0,0,1],'msize',15,'Ax',Ax); %original
Stmp.DrawPointCloud(Stmp.PsuedoInverseCP.CombinePatches(30),'color',[1,1,1],'Ax',Ax); %compact
Stmp.DrawPointCloud([Phd;Qhd],'color',[1,0,0],'msize',20,'Ax',Ax,...
    'title',sprintf('Hausdorff distance %.2g with Radial displacement of %.2g',hd,rhd)); %compact