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
Ax=BezCP.CreateDrawingAxes(fig);

pausetime=0;
% Stmp.RegisteredBezCP.DrawBezierPatches('pausetime',pausetime);
Stmp.StmpBezCP.DrawBezierPatches('Ax',Ax,'pausetime',pausetime);
Stmp.StmpBezCP.DrawControlPoints('Ax',Ax,'pausetime',pausetime);

%% Draw point clouds - orginial, compact and surfaces
Ax=BezCP.CreateDrawingAxes;
H1=BezCP.DrawPointCloud(Stmp.PointCloud,'color',[1,0,0],'msize',20,'Ax',Ax); %original
H2=BezCP.DrawPointCloud(Stmp.Compact,'color',[0,0,1],'msize',20,'Ax',Ax); %compact
H3=BezCP.DrawPointCloud(Stmp.StmpBezCP.Patches2PointCloud(30),'color',[1,1,1],'Ax',Ax); %compact
H4=BezCP.DrawPointCloud(Stmp.StmpBezCP.Vertices,'color',[0,1,0],'Ax',Ax,'msize',20); %CP vertices
lgnd=legend(Ax,'\color{white}Scan','\color{white}compact',...
    '\color{white}bezier surface mesh','\color{white}control points');
set(lgnd,'color',0.2*[1,1,1]);
%% Hausdorff distance with plot
Ax=BezCP.CreateDrawingAxes;
Stmp.HausdorffAsses('Ax',Ax)