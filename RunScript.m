%Calculate 
Stmp=Bez4Stmp('Roie.stl');

%Draw Patches
Stmp.CP.DrawBezierPatches;
Stmp.CP.DrawBezierPatches('curvature','gaussian','PauseTime',0);

%Draw point clouds - orginial, compact and surfaces
Ax=Stmp.DrawPointCloud(Stmp.PointCloud,'color',[1,0,0],'msize',20); %original
Stmp.DrawPointCloud(Stmp.Compact,'color',[0,0,1],'msize',20,'Ax',Ax); %compact
Stmp.DrawPointCloud(Stmp.CP.CombinePatches(30),'color',[1,1,1],'Ax',Ax); %compact
Stmp.DrawPointCloud(Stmp.CP.Vertices,'color',[0,1,0],'Ax',Ax,'msize',20); %CP vertices

%Hausdorff distance with plot
P=Stmp.PointCloud.Location;
Q=Stmp.CP.CombinePatches(30);
szP=size(P); szQ=size(Q);
if numel(szP)==3, P=reshape(P,szP(1)*szP(2),3); end
if numel(szQ)==3, Q=reshape(Q,szQ(1)*szQ(2),3); end
Threshold=30;
P=P(P(:,3)>Threshold,:); Q=Q(Q(:,3)>Threshold,:); %filter buttom noise
[hd,pInd,qInd]=Stmp.Hausdorff(P,Q);
Phd=P(pInd,:); Qhd=Q(qInd,:);
%calculate radial distance (radial off mean point of Phd and Qhd)
m=mean([Phd;Qhd]);
Xcntr=Stmp.Xcenter;
if m(3)>Xcntr(3), t=(m-Xcntr)/norm(m-Xcntr,2);
else, t=[(m(1:2)-Xcntr(1:2)),0]/norm(m(1:2)-Xcntr(1:2),2); end
rhd=dot(Phd-Qhd,t);
%plot
Ax=Stmp.DrawPointCloud(Stmp.PointCloud,'color',[0,0,1],'msize',15); %original
Stmp.DrawPointCloud(Stmp.CP.CombinePatches(30),'color',[1,1,1],'Ax',Ax); %compact
Stmp.DrawPointCloud([Phd;Qhd],'color',[1,0,0],'msize',20,'Ax',Ax,...
    'title',sprintf('Hausdorff distance %.2g with Radial displacement of %.2g',hd,rhd)); %compact