Stmp=Bez4Stmp('Roie.stl');
Ax=Stmp.DrawPointCloud(Stmp.PointCloud);
% Stmp.DrawPointCloud(Stmp.CompactNodes,'Ax',Ax,'color',[0,0,1],'msize',20);
% Stmp.DrawPointCloud(Stmp.BlownNodes,'Ax',Ax,'color',[1,1,1],'msize',20);
Stmp.DrawBezierPatches('Ax',Ax,'PauseTime',0.2);

Stmp.DrawBezierPatches('curvature','gaussian','PauseTime',0);