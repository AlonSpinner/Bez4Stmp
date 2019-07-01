Stmp=Bez4Stmp('Roie.stl');
Ax=Stmp.DrawPointCloud(Stmp.PointCloud);
Stmp.DrawPointCloud(Stmp.CompactControlPoints,'Ax',Ax,'color',[0,0,1],'msize',20);
Stmp.DrawPointCloud(Stmp.BlownControlPoints,'Ax',Ax,'color',[1,1,1],'msize',20);

a=reshape(Stmp.BlownControlPoints.Location,[],3);
b=Stmp.PointCloud.Location;
[hd,aind,bind]=Stmp.Hausdorff(a,b);
Stmp.DrawPointCloud([a(aind,:);b(bind,:)],'Ax',Ax,'color',[1,0,0],'msize',30);