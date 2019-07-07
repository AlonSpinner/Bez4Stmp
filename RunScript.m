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
P=P(P(:,3)>10,:); Q=Q(Q(:,3)>10,:); %filter buttom noise
[hd,pInd,qInd]=Stmp.Hausdorff(P,Q);
Ax=Stmp.DrawPointCloud(Stmp.PointCloud,'color',[0,0,1],'msize',15); %original
Stmp.DrawPointCloud(Stmp.CP.CombinePatches(30),'color',[1,1,1],'Ax',Ax); %compact
Stmp.DrawPointCloud([P(pInd,:);Q(qInd,:)],'color',[1,0,0],'msize',20,'Ax',Ax,...
    'title',sprintf('Hausdorff distance = %.2g',hd)); %compact