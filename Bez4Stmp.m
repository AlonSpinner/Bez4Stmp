%Written by Alon Spinner @ 6-7/19
classdef Bez4Stmp
    properties (SetAccess = private)  %equivalent to "read only" (Getable but not Setable)
        %all properties below are initalized in constructor
        Scan %pointCloud mx3. Point Cloud provided to algorithm
        T2ZScan %pointCloud mx3. Scan after processing Transform2Z
        PointCloud %pointCloud mx3. PointCloud after processing (Transform2Z and downsampling)
        BoundingCylinderLength %%double 1x1.
        BoundingCylinderRadius %%double 1x1.
        SphereRadius %double 1x1.
        Xcenter %double [x,y,z].
        Compact %m x n x 3 double matrix of compact node point cloud
        RegisteredBezCP %BezCP class
        StmpBezCP %BezCP class
    end
    properties (SetAccess = public) %can be changed by user. default values defined here
        GridLengthFcn=@(R,L) R*(0.5*R/L); %after downsampling: VerticesAmount*(GridLength^2)/(2*pi*R*L)~1 by construction
        Cap=true; %true/false (patch on top of sphere or not). if true Slices must be 4.
        SphLayers=2;
        CylLayers=2;
        Slices=4;
        BezierOrder=3;
        XcenterCalculationMethod='normalSTD'; %'GeomtreyAssumption'/'normalSTD'/'normalThreshold'
    end
    methods %public methods, accessible from outside the class
        function obj=Bez4Stmp(PntCldIn,varargin) %Constructor
            %{
            %Input:
            %PntCldIn - points to xyz double matrix of size mx3 point
            %cloud OR pointCloud class object.
            %depending on it's class, assumes the following:
            %double: xyz double matrix of size mx3 [x,y,z]
            %char: file path
            
            %Varargin Input:
            %Time - true/false. measure time of function process (default== true)
            %MaxOptimizationIterions - (default 30)
            
            %Output:
            %initalized class
            %}
            %% Verify and modify inputs
            if nargin==0 %if no input was provided
                clear obj %delete default object which was created when method was called
                error(sprintf(['Input required\n\n',...
                    'Input must be one of the following:\n',...
                    ' -mx3 double matrix containing [x,y,z]\n',...
                    ' -pointCloud class variable\n',...
                    ' -file name of stl\n',...
                    ' -file name of .m file containing pointCloud or mx3 double matrix named "xyz"']));
            end
            
            %obtain xyz - mx3 double class of point cloud vertices
            switch class(PntCldIn)
                case 'double'
                    xyz=PntCldIn;
                case 'char'
                    FileFullName=PntCldIn;
                    [~,FileName,ext]=fileparts(FileFullName);
                    switch ext
                        case '.stl'
                            [~,xyz]=stlread(FileFullName); %obtain vertices of STL
                        case '.mat'
                            s=whos;
                            doubleClassMatches= strcmp({s.class},'double');
                            pointCloudClassMatches= strcmp({s.class},'pointCloud');
                            NameMatches=strcmp({s.name},FileName);
                            if any((doubleClassMatches | pointCloudClassMatches) & NameMatches)
                                eval(['xyz=',s(NameMatches).name,';']);
                            elseif any(pointCloudClassMatches) %take the first pointCloud class variable in list
                                vars={s(pointCloudClassMatches).name};
                                eval(['xyz=',vars{1},';']);
                            elseif any(doubleClassMatches) %take the first double class variable in list
                                vars={s(doubleClassMatches).name};
                                eval(['xyz=',vars{1},';']);
                            end 
                            %convert input to double matrix incase its pointCloud
                            if isa(xyz,'pointCloud'), xyz=xyz.Location; end
                    end
                case 'pointCloud' %incase user went ahead and created a point cloud object
                    xyz=PntCldIn.Location;
                otherwise
                    clear obj %delete default object which was created when method was called
                    error(sprintf(['Wrong Input\n\n',...
                        'Input must be one of the following:\n',...
                        ' -mx3 double matrix containing [x,y,z]\n',...
                        ' -pointCloud class variable\n',...
                        ' -file name of stl\n',...
                        ' -file name of .m file containing pointCloud or mx3 double matrix']));
            end
            
            %default values:
            %rest of default values are defined on public parameters
            Time=true; MaxOptimizationIterions=20; OptDisplay='iter';
            %Obtain inputs
            for ind=1:2:length(varargin)
                comm=lower(varargin{ind});
                switch comm
                    case 'time'
                        Time=varargin{ind+1};
                    case 'maxoptimizationiterions'
                        MaxOptimizationIterions=varargin{ind+1};
                    case 'optdisplay'
                        OptDisplay=varargin{ind+1};
                    case 'gridlengthfcn'
                        obj.GridLengthFcn=varargin{ind+1};
                    case 'cap'
                        obj.Cap=varargin{ind+1};
                    case 'sphlayers'
                        obj.SphLayers=varargin{ind+1};
                    case 'cyllayers'
                        obj.CylLayers=varargin{ind+1};
                    case 'slices'
                        obj.Slices=varargin{ind+1};
                    case 'bezierorder'
                        obj.BezierOrder=varargin{ind+1};
                    case 'XcenterCalculationMethod'
                        obj.XcenterCalculationMethod=varargin{ind+1};
                end
            end
            
            %% Run
            if Time, tic; end %start counting time of process
            
            %introduce to parameters before any preprocessing to object
            obj.Scan=pointCloud(xyz);
            %align with z axis, sets base on z=0, sort vertices by z value
            T2Z_xyz=Transform2Z(xyz);
            %find Radius and Length of bounding Cylinder
            [R,L]=BoundingCylinder(T2Z_xyz);
            obj.BoundingCylinderLength=L; %introduce to parameters
            obj.BoundingCylinderRadius=R; %introduce to parameters
            
            %Create pointCloud object and introduce to parameters
            PntCld=pointCloud(T2Z_xyz); %create point cloud object from preprocessed point cloud
            obj.T2ZScan=PntCld; %introduce to parameters
            
            %DownSample - final part of preprocessing - and introduce to
            %parameters
            GridLength=obj.GridLengthFcn(R,L); %obtain voxel edge length by default function
            PntCld=pcdownsample(PntCld,'gridAverage',GridLength); %Downsample by grid legnth
            obj.PointCloud=PntCld; %introduce to parameters
            
            %Find Xcenter - point seperating between sphereical and
            %Cylinderical part
            Xcntr=FindXcenter(PntCld,obj.XcenterCalculationMethod);
            obj.Xcenter=Xcntr; %introduce to parameters
            
            %Find radius of sphereical part
            r=L-Xcntr(:,3);
            obj.SphereRadius=r;
            
            %Create Compact Control Points nodes
            minR=min(r,R); maxR=max(r,R);
            [Ncyl,Nsph,Ncircum]=MeshNodesAmount(obj.SphLayers,obj.CylLayers,obj.Slices,obj.BezierOrder);
            P=CircleBellCrv(minR,L,Ncyl,Nsph,'RemoveTopEdge',obj.Cap);
            P=flipud(P); %we want highest z in the first row for BezCP construction
            xyz=RotateXZcurve(P,Ncircum); %rotate to create double size Nvert x Ncircum x 3
            xyz(:,:,1)=xyz(:,:,1)+Xcntr(1); xyz(:,:,2)=xyz(:,:,2)+Xcntr(2); %align vertical axis to pass through Xcenter
            obj.Compact=xyz;
            
            %Create compact BezCP (CompactCP)
            if obj.Cap, CompactCP=BezCP(xyz,obj.BezierOrder,'method','CircularWithCap');
            else,  CompactCP=BezCP(xyz,obj.BezierOrder,'method','Circular'); end
            
            %Optimize CompactCP by radially expanding the vertices to PntCld
            h=helpdlg('Optimization occuring. Please wait. . .');
            obj.RegisteredBezCP=RadialOptimization4CP(CompactCP,PntCld,[-minR,+maxR],Xcntr,'Display',OptDisplay,...
                'MaxIterations',MaxOptimizationIterions); %expand CompactCP to PntCld
            close(h); %close helpdlg
            
            %find control points describiing optimized CP.Vertices mesh surface
            obj.StmpBezCP=obj.RegisteredBezCP.PesudoInverseVertices;
            
            %if no cap, top edges whose location may vary across
            %top patches due to optimization to produce a peak point
            if ~obj.Cap
                obj.StmpBezCP=obj.StmpBezCP.UnifyEdge2Peak('top');
            end
            
            if Time, toc; end %display time it took for process
        end
        function obj=UpdateObj(obj,varargin)
            %{
            %Varargin Input:
            %Time - true/false. measure time of function process (default== true)
            %MaxOptimizationIterions - default (30)
            
            %Output:
            %obj - Updates the private access properties after user had changed
            %the public access properties
                %Private access properties updated:
                    %PointCloud  (reminder: GridLengthFcn changes)
                    %Xcenter     (reminder: Xcenter calculation method can)
                    %change
                    %SphereRadius
                    %Compact
                    %CP
            %}
            
            %default values:
            Time=true; MaxOptimizationIterions=20; OptDisplay='iter';
            %Obtain inputs
            for ind=1:2:length(varargin)
                comm=lower(varargin{ind});
                switch comm
                    case 'time'
                        Time=varargin{ind+1};
                    case 'maxoptimizationiterions'
                        MaxOptimizationIterions=varargin{ind+1};
                    case 'display'
                        OptDisplay=varargin{ind+1};
                end
            end
            %% Run
            if Time, tic; end %start counting time of process
            
            %introduce to parameters before any preprocessing to object
            obj.Scan=pointCloud(xyz);
            %align with z axis, sets base on z=0, sort vertices by z value
            T2Z_xyz=Transform2Z(xyz);
            %find Radius and Length of bounding Cylinder
            [R,L]=BoundingCylinder(T2Z_xyz);
            obj.BoundingCylinderLength=L; %introduce to parameters
            obj.BoundingCylinderRadius=R; %introduce to parameters
            
            %Create pointCloud object and introduce to parameters
            PntCld=pointCloud(T2Z_xyz); %create point cloud object from preprocessed point cloud
            obj.T2ZScan=PntCld; %introduce to parameters
            
            %DownSample - final part of preprocessing - and introduce to
            %parameters
            GridLength=obj.GridLengthFcn(R,L); %obtain voxel edge length by default function
            PntCld=pcdownsample(PntCld,'gridAverage',GridLength); %Downsample by grid legnth
            obj.PointCloud=PntCld; %introduce to parameters
            
            %Find Xcenter - point seperating between sphereical and
            %Cylinderical part
            Xcntr=FindXcenter(PntCld,obj.XcenterCalculationMethod);
            obj.Xcenter=Xcntr; %introduce to parameters
            
            %Find radius of sphereical part
            r=L-Xcntr(:,3);
            obj.SphereRadius=r;
            
            %Create Compact Control Points nodes
            minR=min(r,R); maxR=max(r,R);
            [Ncyl,Nsph,Ncircum]=MeshNodesAmount(obj.SphLayers,obj.CylLayers,obj.Slices,obj.BezierOrder);
            P=CircleBellCrv(minR,L,Ncyl,Nsph,'RemoveTopEdge',obj.Cap);
            P=flipud(P); %we want highest z in the first row for BezCP construction
            xyz=RotateXZcurve(P,Ncircum); %rotate to create double size Nvert x Ncircum x 3
            xyz(:,:,1)=xyz(:,:,1)+Xcntr(1); xyz(:,:,2)=xyz(:,:,2)+Xcntr(2); %align vertical axis to pass through Xcenter
            obj.Compact=xyz;
            
            %Create compact BezCP (CompactCP)
            if obj.Cap, CompactCP=BezCP(xyz,obj.BezierOrder,'method','CircularWithCap');
            else,  CompactCP=BezCP(xyz,obj.BezierOrder,'method','Circular'); end
            
            %Optimize CompactCP by radially expanding the vertices to PntCld
            h=helpdlg('Optimization occuring. Please wait. . .');
            obj.RegisteredBezCP=RadialOptimization4CP(CompactCP,PntCld,[-minR,+maxR],Xcntr,'Display',OptDisplay,...
                'MaxIterations',MaxOptimizationIterions); %expand CompactCP to PntCld
            close(h); %close helpdlg
            
            %find control points describiing optimized CP.Vertices mesh surface
            obj.StmpBezCP=obj.RegisteredBezCP.PesudoInverseVertices;
            
            %if no cap, top edges whose location may vary across
            %top patches due to optimization to produce a peak point
            if ~obj.Cap
                obj.StmpBezCP=obj.StmpBezCP.UnifyEdge2Peak('top');
            end
            
            if Time, toc; end %display time it took for process
        end
        function varargout=HausdorffAsses(obj,varargin)
            %{
            Varagin inputs:
            Ax - handle of axes. if provided - drawns on them. if not -
            doesnt draw.
            N - number of points drawn in each evauluated patch are N^2. by default N=30
            ZThreshold - filters out all points below ZThreshold
            
            Outputs by varargout:
            varargout{1} = hausdorff distance
            varargout{2} = radial hausdorff distance
            varargout{3} = handle array of graphic objects
            %}
            N=30; zThreshold=30;
            for ind=1:2:length(varargin)
                comm=lower(varargin{ind});
                switch comm
                    case 'n'
                        N=varargin{ind+1};
                    case 'ax'
                        Ax=varargin{ind+1};
                    case 'zthreshold'
                        zThreshold=varargin{ind+1};
                end
            end
            
            %find hausdroff distance and the points who make it
            P=obj.PointCloud.Location;
            P=P(P(:,3)>zThreshold,:); %filter out buttom noise
            Q=obj.StmpBezCP.CombinePatches(N);
            szQ=size(Q);
            Q=reshape(Q,szQ(1)*szQ(2),3); %reshape Q to mx3
            Q=Q(Q(:,3)>zThreshold,:); %filter buttom noise
            [hd,pInd,qInd]=BezCP.Hausdorff(P,Q); %computing hausdorff distance
            Phd=P(pInd,:); Qhd=Q(qInd,:); %obtain points whose distance is the hausdorff distance
            
            %calculate radial distance (radial off mean point of Phd and Qhd)
            m=mean([Phd;Qhd]);
            Xcntr=obj.Xcenter;
            if m(3)>Xcntr(3), t=(m-Xcntr)/norm(m-Xcntr,2);
            else, t=[(m(1:2)-Xcntr(1:2)),0]/norm(m(1:2)-Xcntr(1:2),2); end
            rhd=abs(dot(Phd-Qhd,t)); %radial hausdorff distance
            
            %plot if axes were provided
            if exist('Ax','var') && isgraphics(Ax,'Axes')
                Handles=gobjects(1,3);
                Handles(1)=BezCP.DrawPointCloud(obj.PointCloud,'color',[0,0,1],'msize',15,'Ax',Ax); %original
                Handles(2)=BezCP.DrawPointCloud(obj.StmpBezCP.CombinePatches(30),'color',[1,1,1],'Ax',Ax); %compact
                Handles(3)=BezCP.DrawPointCloud([Phd;Qhd],'color',[1,0,0],'msize',20,'Ax',Ax,...
                    'title',sprintf('Hausdorff distance %.2g with Radial displacement of %.2g',hd,rhd)); %compact
            end
            
            %Create output
            if nargout==0, varargout={};
            else, varargout={hd,rhd,Handles}; end
        end
    end
end
%% algorithm pipeline functions
function Newxyz=Transform2Z(xyz)
%{
Translate xyz points so that CGxy=[0,0], floor  is z=0, and primary
axis of the convex shape will be the Z axis.
NewXYZ is also sorted by Z (first row with minimal Z and last row with
maximal)

Input:
xyz - [x,y,z] numeric matrix mX3

Output:
Newxyz - [x,y,z] numeric matrix mX3.

centerlize xyz and compute pca analysis
%}
XYZ0=xyz-mean(xyz);
V=pca(XYZ0);
t=V(:,1); %col vec

%define rotation matrix that rotates t->z by rodriguez formula
% https://math.stackexchange.com/questions/180418/calculate-rotation-matrix-to-align-vector-a-to-vector-b-in-3d
% Rodrigues's rotation formula gives the result of a rotation of a vector a around a rotation axis k by the angle ?.
% We can make use of this by realizing that, to rotate a unit vector t into z, we simply need to rotate a by ? around (t+z)/2.
% With this, one gets the beautiful
z=[0;0;1]; %col vec
Rtz=2*((t+z)*(t+z)')/((t+z)'*(t+z))-eye(3);
R=zeros(4,4); R(1:3,1:3)=Rtz; R(4,4)=1;

%Rotate the data so the main axis is parallel to z
OnesVec=ones(size(XYZ0,1),1);
Rxyz=(R*[XYZ0,OnesVec]')';
Rxyz=Rxyz(:,1:3);

%translate the data so we start from [0,0,0]
Txyz=Rxyz(:,1:3); %get rid of the ones vector
Txyz(:,3)=Txyz(:,3)-min(Rxyz(:,3)); %translate

%Sort by Z
[~,Indx]=sort(Txyz(:,3));
Newxyz=Txyz(Indx,:);
end
function [R,L]=BoundingCylinder(xyz)
%{
Assumption:
xyz was already processed by function Transfom2Z:
it is aligned with z axis, has base on z=0

Input
xyz - matrix of [x,y,z] point cloud ((mx3)

Output:
R - radius of bounding cylinder
L - Length of bouding cylinder

Calculate L
%}

L=max(xyz(:,3));

% Calculate R
% algerbric error minimization (not geometric fit)
%
% Reference: "Least-squares fitting of circles and ellipses", W. Gander,
% G. Golub, R. Strebel - BIT Numerical Mathematics, 1994, Springer

% This implementation copyright Richard Brown, 2007, but is freely
% available to copy, use, or modify as long as this line is maintained

%project all XYZ points to XY plane and find convex hull.
%Circl will be  fitted on convex hull
xy=xyz(:,[1,2]);
CHull=xy(convhull(xy),:);  %[x,y] mat. find convex hull points

[~,R]=fitcircle(CHull,'linear');
end
function Xcenter=FindXcenter(PtCld,method,varargin)
%{
Assumption:
--PtCld was already processed by function Transfom2Z:
it is aligned with z axis, has base on z=0
--PtCld has already gone through downsampling using GridLengthFcn(R,L) and
function pcdownsample.

Input:
PtCldP  - pointCloud object
Method - 'GeomtreyAssumption'/'normalSTD'/'normalThreshold'
Varargin - depending on method:
normalSTD - BoundingCylinderLength
GeomtreyAssumption: BoundingCylinderRadius
normalThreshold: requires Threshold scalar ranges from [0,1]
suggestion: 0.4

Output:
Xcenter [x,y,z] is the node that seprates between the "sphere" part and the
"cylinder" part in the given cloud point

Assumptions:
CldP main axis
    ___
  /     \
/        |
|     x   |
|         |
|         |
|         |
%}

for ind=1:2:length(varargin)
    comm=lower(varargin{ind});
    switch comm
        case 'boundingcylinderradius'
            R=varargin{ind+1};
        case 'Threshold'
            Threshold=varargin{ind+1};
    end
end

switch method
    case 'GeomtreyAssumption'
        if ~exist('R','var') || isempty(R)
            error('BoundingCylinderRadius not introduced to function');
        end
        Znode=L-R;
    case 'normalSTD'
        %checks for biggest change in standart variation in transient
        %signal and determines it as the point of interest.
        
        %we check normals with 8 nearest neighbors which is more or less
        %like checking for vertices within sphere with radius GridLength.
        %*|*|*
        %*|#|*
        %*|*|*
        
        %sometimes there is noise in the lower parts of the given scan that prodces vertical normals
        zmax=max(PtCld.Location(:,3));
        pass=find(PtCld.Location(:,3)>zmax/4); %all vertices that pass the condition
        CleanPtCld=select(PtCld,pass);
        k=8;
        normals=pcnormals(CleanPtCld,k);
        Idx=findchangepts(normals(:,3),'Statistic','std'); %check for the point in which the std changes the moves
        Znode=CleanPtCld.Location(Idx,3);
    case 'normalThreshold'
        if ~exist('Threshold','var') || isempty(Threshold) || (Threshold>1 || Threshold<0)
            error('Failed to compute Xcenter. Need to provide a Threshold in range [0,1] with the chosen method')
        end
        %sometimes there is noise in the lower parts of the given scan that prodces vertical normals
        zmax=max(PtCld.Location(:,3));
        pass=find(PtCld.Location(:,3)>zmax/4); %all vertices that pass the condition
        CleanPtCld=select(PtCld,pass);
        k=8; %see explination in normalSTD method.
        normals=pcnormals(CleanPtCld,k);
        pass=find(abs(normals)>Threshold); %all vertices that pass the condition
        Idx=pass(1); %take first vertex that pass condition
        Znode=CleanPtCld.Location(Idx,3);
end

%find [x,y] of Xcenter by mean on all the points above it (higher than Znode)
pass=find(PtCld.Location(:,3)>Znode); %all vertices that pass the condition
AbvPtCld=select(PtCld,pass);
Xcenter=[mean(AbvPtCld.Location(:,[1,2])),Znode];

end
function [Ncyl,Nsph,Ncircum]=MeshNodesAmount(SphLayers,CylLayers,Slices,BezierOrder)
%{
Input:
CylLayers, SphLayers - amount of layers in spherical/cylinderical part
Slices - amount of bezier patches on the circumference in a given
Layer
BezierOrder - polynomial order of bezier patches (symmeterical for both
parameters)

Output:
Ncyl - amount of points in a vertical line (defacto z direction) in
cylinderical part
Nsph - amount of points in a vertical line (defacto z direction) in
spherical part
Ncircm - amount of points in a horizontal (circumference of stump) line

Ncircum=(BezierOrder+1)+(Slices-1)*BezierOrder=CakeLayers*BezierOrder
Nvert=(BezierOrder+1)+(Layers-2)*BezierOrder+(BezierOrder-1)=CakeSlices*BezierOrder+1
%}

Ncyl=CylLayers*BezierOrder+1;
Nsph=SphLayers*BezierOrder; %1 less point than Ncyl as it is shared with cylinder
Ncircum=Slices*BezierOrder;
end
function P=CircleBellCrv(R,h,Ncyl,Nsph,varargin)
%{
Input:
R - radius
h - height
Ncyl - number of points in "cylinderical" part
Nsph - number of points in "Spherical" part

Varargin Inputs:
RemoveTopEdge - true/false of spherical part. (default=true)
false ensures that the spherical top points are all the
same point (the peak point is repeated)
RemoveButtomEdge - true/false of spherical part (default=true)
false ensures that spherical buttom points and cylinderical
points are the same (repeated points)
Output:
P - [x,y] (size (Nsph+Ncyl)x2 of evaluated points
%}

%---           ^
%  ^  \        |
%  |    \      |
%  R      \    |
%  |       \   |
%  v       |   h
%          |   |
%          |   |
% <---R--->|   v


%default values:
RemoveButtomEdge=true; RemoveTopEdge=true;
%Obtain inputs
for ind=1:2:length(varargin)
    comm=lower(varargin{ind});
    switch comm
        case 'removebuttomedge'
            RemoveButtomEdge=varargin{ind+1};
        case 'removetopedge'
            RemoveTopEdge=varargin{ind+1};
    end
end

%built function. note:
%((z-(h-R))/R)=sin(theta)--->sqrt(1-((z-(h-R))/R)^2)=cos(theta)
fcnrz=@(z) R-heaviside(z-(h-R))*R*(1-sqrt(1-((z-(h-R))/R)^2));

%build vectors
zcyl=linspace(0,h-R,Ncyl);
zsph=linspace(h-R,h,Nsph+RemoveTopEdge+RemoveButtomEdge);
zsph=zsph(1+RemoveButtomEdge:end-RemoveTopEdge);
z=[zcyl,zsph];

%evalulate r as a function of z
r=arrayfun(fcnrz, z);

%build output
P=[r(:),z(:)];
end
function xyz=RotateXZcurve(xz,n)
%{
Input:
xz [x,z] of evalualted points on 2D curve (mx2 array)
n - multiplyer amount

Output:
xyz  evaluated points on 3D surface (mxn3 array) where depth dimension is
      orginized [x,y]z.
      created by turnning the xy curve n times around the y axis. each
      time by 2*pi/N amount
%}

m=size(xz,1); %amount of evaluated points on curve
theta=linspace(0,2*pi,n+1);
theta=theta(2:end); %first and last point were the same
x=xz(:,1); z=xz(:,2);
[X,Y,Z]=deal(zeros(m,n)); %initalize
for i=1:n
    X(:,i)=x*cos(theta(i));
    Y(:,i)=x*sin(theta(i));
    Z(:,i)=z;
end
xyz=cat(3,X,Y,Z);
end
function CP=RadialOptimization4CP(CP,StaticPntCld,Rbounds,Xcenter,varargin)
%{
Blow up CP.Vertices to StaticPntCld using fmincon (optimization)

two methods for radial optimization in this code:
--cylinderical: radial optimization from the Xcenter axis from points
below Xcenter (zwise)
--spherical: radial optimization from Xcenter for points above Xcenter
(zwise)

Algorithm Assumption:
--CP.Vertices are bounded within StaticPntCld.
--to register bezier surfaces onto StaticPntCLd one must expand/shrink CP.Vertices radially
from a given point/axis.

Example for expansion from a given axis:
|  <-|   |->  |
|  <-|   |->  |
|  <-|   |->  |

Input:
StaticPntCld - Point cloud in format of [x,y,z] (mx3) OR pointCloud
object with similar .Location value
CP.v - see notes below
OR pointCloud object with similar.Location value
RBounds - [Minimal radius,Maximal radius] to add to a point in MovPntCld
Xcenter - [x,y,z] mx3 point in 3D space. for spherical mehthod it
is the point from which radial expansion happens about. for cylinderical
method Xcenter(1,2) is the [x,y] coordante where the centeral axis passes
through and expansion happens about

Varargin inputs:
DiffMinChange/DiffMaxChange - double, scalar. min/max change in radial
values of points
MaxIterations - integer,scalar.
Display - 'none'/'iter/'final'. displays data in fmincon. see fmincon
options for more

Output:
CP - updated, optimized CP
%}

%varargin input and their defaults
DeltaR=Rbounds(2)-Rbounds(1);
DiffMinChange=DeltaR/1000; DiffMaxChange=DeltaR/10; MaxIterations=20;
Display='none';
for ind=1:2:length(varargin)
    comm=lower(varargin{ind});
    switch comm
        case 'diffminchange'
            DiffMinChange=varargin{ind+1};
        case 'diffmaxchange'
            DiffMaxChange=varargin{ind+1};
        case 'maxiterations'
            MaxIterations=varargin{ind+1};
        case 'display'
            Display=varargin{ind+1};
    end
end

%if point clouds give in uncomfi formats, change them
%StaticPntCld needs to be pointCloud object for cost function (use of
%findNearestNeighbors)
if ~isa(StaticPntCld,'pointCloud'), StaticPntCld=pointCloud(StaticPntCld); end

V=CP.Vertices; %to shorten writing
m=size(V,1); %find amount of points to move

%find t - matrix mx3 of radial direction per point
Isph=V(:,3)>Xcenter(3); %boolean array pointing to spherical part
Icyl=~Isph; %boolean array pointing to cylinderical part
tsph=(V-Xcenter)./vecnorm(V-Xcenter,2,2); %sphereical radial vectors
tcyl=[V(:,[1,2])-Xcenter(1:2),zeros(m,1)]./vecnorm(V(:,[1,2])-Xcenter(1:2)+eps,2,2); %cylinderical radial vectors
t=Isph.*tsph+Icyl.*tcyl; %radial vectors of all CP.Vertices depending on their classification
%note:
%tcyl denominator has +eps in it to avoid nan values as
%vecnorm(V(:,[1,2])-Xcenter(1:2),2,2) may result in zero

%Optimization
fun=@(r) CostFcn(CP,StaticPntCld,t,r);
options = optimoptions('fmincon','Algorithm','interior-point','Display',Display,...
    'DiffMaxChange',DiffMaxChange,'DiffMinChange',DiffMinChange,...
    'MaxIterations',MaxIterations,'MaxFunctionEvaluations',MaxIterations*m);
lb=Rbounds(1)*ones(m,1);
ub=Rbounds(2)*ones(m,1);
r0=zeros(m,1);
rOpt=fmincon(fun,r0,[],[],[],[],lb,ub,[],options);

%Create output
CP.Vertices=V+t.*rOpt;
    function Cost=CostFcn(CP,StaticPntCld,t,r)
        %Input:
        %CP - BezCP class
        %StaticPntCld  - pointCloud class
        %t - vectors size mx3 of radial direction.
        %r - vector mx1 of radius values (optimization parameter)
        
        %Output:
        %Cost=least squares cost essentially (CP.Vertices+t.*r - StaticPntCld)^2
        
        %note: Nested functions variables who share a common name with
        %their parent function will become global. for such reason, 'N' was
        %used instead of 'm'
        Cost=0;
        CP.Vertices=CP.Vertices+t.*r;
        for i=1:size(CP.Vertices,1)
            [~,d]=findNearestNeighbors(StaticPntCld,CP.Vertices(i,:),1);
            Cost=Cost+d^2;
        end
    end
end
%% extra functions (note: can't call one method from another)
%read stl
function varargout=stlread(file)
%{
Works only with binary STLS
imported from https://www.mathworks.com/matlabcentral/fileexchange/22409-stl-file-reader
slightly edited to remove ASCII parts which did not work
anyhow

STLREAD imports geometry from an STL file into MATLAB.
   FV = STLREAD(FILENAME) imports triangular faces from the ASCII or binary
   STL file idicated by FILENAME, and returns the patch struct FV, with fields
   'faces' and 'vertices'.

   [F,V] = STLREAD(FILENAME) returns the faces F and vertices V separately.

   [F,V,N] = STLREAD(FILENAME) also returns the face normal vectors.

   The faces and vertices are arranged in the format used by the PATCH plot
   object.
Copyright 2011 The MathWorks, Inc.
%}
if ~exist(file,'file')
    error(['File ''%s'' not found. If the file is not on MATLAB''s path' ...
        ', be sure to specify the full path to the file.'], file);
end

fid = fopen(file,'r');
if ~isempty(ferror(fid))
    error(lasterror); %#ok
end

M = fread(fid,inf,'uint8=>uint8');
fclose(fid);

[f,v,n] = stlbinary(M);

varargout = cell(1,nargout);
switch nargout
    case 2
        varargout{1} = f;
        varargout{2} = v;
    case 3
        varargout{1} = f;
        varargout{2} = v;
        varargout{3} = n;
    otherwise
        varargout{1} = struct('faces',f,'vertices',v);
end
end
function [F,V,N]=stlbinary(M)
F = [];
V = [];
N = [];

if length(M) < 84
    error('MATLAB:stlread:incorrectFormat', ...
        'Incomplete header information in binary STL file.');
end

% Bytes 81-84 are an unsigned 32-bit integer specifying the number of faces
% that follow.
numFaces = typecast(M(81:84),'uint32');
%numFaces = double(numFaces);
if numFaces == 0
    warning('MATLAB:stlread:nodata','No data in STL file.');
    return
end

T = M(85:end);
F = NaN(numFaces,3);
V = NaN(3*numFaces,3);
N = NaN(numFaces,3);

numRead = 0;
while numRead < numFaces
    % Each facet is 50 bytes
    %  - Three single precision values specifying the face normal vector
    %  - Three single precision values specifying the first vertex (XYZ)
    %  - Three single precision values specifying the second vertex (XYZ)
    %  - Three single precision values specifying the third vertex (XYZ)
    %  - Two unused bytes
    i1    = 50 * numRead + 1;
    i2    = i1 + 50 - 1;
    facet = T(i1:i2)';
    
    n  = typecast(facet(1:12),'single');
    v1 = typecast(facet(13:24),'single');
    v2 = typecast(facet(25:36),'single');
    v3 = typecast(facet(37:48),'single');
    
    n = double(n);
    v = double([v1; v2; v3]);
    
    % Figure out where to fit these new vertices, and the face, in the
    % larger F and V collections.
    fInd  = numRead + 1;
    vInd1 = 3 * (fInd - 1) + 1;
    vInd2 = vInd1 + 3 - 1;
    
    V(vInd1:vInd2,:) = v;
    F(fInd,:)        = vInd1:vInd2;
    N(fInd,:)        = n;
    
    numRead = numRead + 1;
end
end %for stlread
%fit circle to find BoundingCylinder R
function [z,r,residual]=fitcircle(x,varargin)
%{
Obtained from "fitcircle.m" from mathworks file exchange
https://www.mathworks.com/matlabcentral/fileexchange/15060-fitcircle-m

FITCIRCLE    least squares circle fit

  [Z, R] = FITCIRCLE(X) fits a circle to the N points in X minimising
  geometric error (sum of squared distances from the points to the fitted
  circle) using nonlinear least squares (Gauss Newton)
      Input
          X : 2xN array of N 2D points, with N >= 3
      Output
          Z : center of the fitted circle
          R : radius of the fitted circle

  [Z, R] = FITCIRCLE(X, 'linear') fits a circle using linear least
  squares minimising the algebraic error (residual from fitting system
  of the form ax'x + b'x + c = 0)

  [Z, R] = FITCIRCLE(X, Property, Value, ...) allows parameters to be
  passed to the internal Gauss Newton method. Property names can be
  supplied as any unambiguous contraction of the property name and are
  case insensitive, e.g. FITCIRCLE(X, 't', 1e-4) is equivalent to
  FITCIRCLE(X, 'tol', 1e-4). Valid properties are:

      Property:                 Value:
      --------------------------------
      maxits                    positive integer, default 100
          Sets the maximum number of iterations of the Gauss Newton
          method

      tol                       positive constant, default 1e-5
          Gauss Newton converges when the relative change in the solution
          is less than tol

  [X, R, RES] = fitcircle(...) returns the 2 norm of the residual from
  the least squares fit

  Example:
      x = [1 2 5 7 9 3; 7 6 8 7 5 7];
      % Get linear least squares fit
      [zl, rl] = fitcircle(x, 'linear')
      % Get true best fit
      [z, r] = fitcircle(x)

  Reference: "Least-squares fitting of circles and ellipses", W. Gander,
  G. Golub, R. Strebel - BIT Numerical Mathematics, 1994, Springer

This implementation copyright Richard Brown, 2007, but is freely
available to copy, use, or modify as long as this line is maintained
%}

error(nargchk(1, 5, nargin, 'struct'))

% Default parameters for Gauss Newton minimisation
params.maxits = 100;
params.tol    = 1e-5;

% Check x and get user supplied parameters
[x, fNonlinear, params] = parseinputs(x, params, varargin{:});

% Convenience variables
m  = size(x, 2);
x1 = x(1, :)';
x2 = x(2, :)';


% 1) Compute best fit w.r.t. algebraic error using linear least squares
%
% Circle is represented as a matrix quadratic form
%     ax'x + b'x + c = 0
% Linear least squares estimate found by minimising Bu = 0 s.t. norm(u) = 1
%     where u = [a; b; c]

% Form the coefficient matrix
B = [x1.^2 + x2.^2, x1, x2, ones(m, 1)];

% Least squares estimate is right singular vector corresp. to smallest
% singular value of B
[U, S, V] = svd(B);
u = V(:, 4);

% For clarity, set the quadratic form variables
a = u(1);
b = u(2:3);
c = u(4);

% Convert to centre/radius
z = -b / (2*a);
r = sqrt((norm(b)/(2*a))^2 - c/a);

% 2) Nonlinear refinement to miminise geometric error, and compute residual
if fNonlinear
    [z, r, residual] = fitcircle_geometric(x, z, r);
else
    residual = norm(B * u);
end

% END MAIN FUNCTION BODY

% NESTED FUNCTIONS
    function [z, r, residual] = fitcircle_geometric(x, z0, r0)
        % Use a simple Gauss Newton method to minimize the geometric error
        fConverged = false;
        
        % Set initial u
        u     = [z0; r0];
        
        % Delta is the norm of current step, scaled by the norm of u
        delta = inf;
        nIts  = 0;
        
        for nIts = 1:params.maxits
            % Find the function and Jacobian
            [f, J] = sys(u);
            
            % Solve for the step and update u
            h = -J \ f;
            u = u + h;
            
            % Check for convergence
            delta = norm(h, inf) / norm(u, inf);
            if delta < params.tol
                fConverged = true;
                break
            end
        end
        
        if ~fConverged
            warning('fitcircle:FailureToConverge', ...
                'Gauss Newton iteration failed to converge');
        end
        z = u(1:2);
        r = u(3);
        f = sys(u);
        residual = norm(f);
        
        
        function [f, J] = sys(u)
            %SYS   Nonlinear system to be minimised - the objective
            %function is the distance to each point from the fitted circle
            %contained in u
            
            % Objective function
            f = (sqrt(sum((repmat(u(1:2), 1, m) - x).^2)) - u(3))';
            
            % Jacobian
            denom = sqrt( (u(1) - x1).^2 + (u(2) - x2).^2 );
            J = [(u(1) - x1) ./ denom, (u(2) - x2) ./ denom, repmat(-1, m, 1)];
        end % sys
        
    end % fitcircle_geometric

% END NESTED FUNCTIONS

end
function [x,fNonlinear,params]=parseinputs(x,params,varargin)
% Make sure x is 2xN where N > 3
if size(x, 2) == 2
    x = x';
end

if size(x, 1) ~= 2
    error('fitcircle:InvalidDimension', ...
        'Input matrix must be two dimensional')
end

if size(x, 2) < 3
    error('fitcircle:InsufficientPoints', ...
        'At least 3 points required to compute fit')
end

% determine whether we are measuring geometric error (nonlinear), or
% algebraic error (linear)
fNonlinear = true;
switch length(varargin)
    % No arguments means a nonlinear least squares with defaul parameters
    case 0
        return
        
        % One argument can only be 'linear', specifying linear least squares
    case 1
        if strncmpi(varargin{1}, 'linear', length(varargin{1}))
            fNonlinear = false;
            return
        else
            error('fitcircle:UnknownOption', 'Unknown Option')
        end
        
        % Otherwise we're left with user supplied parameters for Gauss Newton
    otherwise
        if rem(length(varargin), 2) ~= 0
            error('fitcircle:propertyValueNotPair', ...
                'Additional arguments must take the form of Property/Value pairs');
        end
        
        % Cell array of valid property names
        properties = {'maxits', 'tol'};
        
        while length(varargin) ~= 0
            property = varargin{1};
            value    = varargin{2};
            
            % If the property has been supplied in a shortened form, lengthen it
            iProperty = find(strncmpi(property, properties, length(property)));
            if isempty(iProperty)
                error('fitcircle:UnkownProperty', 'Unknown Property');
            elseif length(iProperty) > 1
                error('fitcircle:AmbiguousProperty', ...
                    'Supplied shortened property name is ambiguous');
            end
            
            % Expand property to its full name
            property = properties{iProperty};
            
            switch property
                case 'maxits'
                    if value <= 0
                        error('fitcircle:InvalidMaxits', ...
                            'maxits must be an integer greater than 0')
                    end
                    params.maxits = value;
                case 'tol'
                    if value <= 0
                        error('fitcircle:InvalidTol', ...
                            'tol must be a positive real number')
                    end
                    params.tol = value;
            end % switch property
            varargin(1:2) = [];
        end % while
        
end % switch length(varargin)

end %parseinputs %for fitcircle
%% functions not in use - ideas not implemented
function [P,Q]=BezBellCurve(R,d,L,Nq)
%{
%Builds a 2D bezier curve and evalutes Nq points on it
%(eqidistance on parameter q in [0,1])

%    nodes (symbolised by "@") are set as such:
%  @4<----R--->@3
%              ^      ^
%              |      |
%              | d    |
%              |      |
% y            v      |
% ^            @2     |
% |            |      |
% |            |      |
% |            |L-d   |  L
%              |      |
%              |      |
%              |      |
%              v      |
%              @1      v   ---> x

%Input:
%R - radius
%d - legth
%L - length
%Nq - number of evaluated points on curve

%Output:
%P - evaluated points on Bezier Curve [x,y,z] (mx2)
%Q - nodes [x,y,z] (mx2)

%Example:
% L=10; R=3; d=3; Nq=10;
% [P,Q]=Bez4Stmp.BellCurve(R,d,L,Nq);
% x=axes; hold(x,'on'); grid(x,'on'); axis(x,'equal');
% scatter(P(:,1),P(:,2),5,'b');
% scatter(Q(:,1),Q(:,2),10,'r');
% %NOTE: d=R, d=0 and d=L are the three interesting cases
%}

%build nodes and parameter vector
Q=[R,0;
    R,L-d;
    R,L;
    0,L];
q=linspace(0,1,Nq);

%Evaluate points
P=zeros(Nq,2);
for i=1:Nq
    P(i,:)=EvalBezCrv_DeCasteljau(Q,q(i));
end

    function R=EvalBezCrv_DeCasteljau(Q,q)
        %Evaluates Bezier Curve by given nodes and parameter
        %value
        
        %Q - nodes in format [x] or [x,y,z] dimension matrix. top row is q=0.
        %q - running parameter of bezier curve. 0=<q<=1
        %R - [x] or [x,y,z] format. point on bezier curve
        % https://pages.mtu.edu/~shene/COURSES/cs3621/NOTES/spline/Bezier/de-casteljau.html
        
        n=size(Q,1)-1; %degree of bezier polynomial
        for k=1:(n+1)
            for j=1:(n+1-k) %doesnt enter in first iteration
                Q(j,:)=(1-q)*Q(j,:)+q*Q(j+1,:);
            end
        end
        R=Q(1,:);
    end

end
function [xyz,Psi,Theta]=CalcSphere(R,Npsi,Ntheta,varargin)
%{
%INPUT:
%R - raidus
%Ntheta - amount of theta points (descrization). theta ranges [0,2pi]
%Npsi - amount of psi points (descrization). psi ranges [0,pi]

%Varagin:
%Trans - translation [x,y,z]
%Half - true/false. true - plots half dome (pointing up)
%Edges - true/false. true - divide psi differently in range [0,pi] so it wont have
%the values of 0,pi

%OUTPUT:
%xyz - size Npsi x Ntheta x 3 where the depth dimension is constructed as [x,y,z]
%Theta, Psi. matrices of size Ntheta x Npsi
%}

%default values:
Half=false; Edges=true;
%Obtain inputs
for ind=1:2:length(varargin)
    comm=lower(varargin{ind});
    switch comm
        case 'trans'
            Trans=varargin{ind+1};
        case 'half'
            Half=varargin{ind+1};
        case 'edges'
            Edges=varargin{ind+1};
    end
end

%create parameter vectors
theta=linspace(0,2*pi,Ntheta+1);
theta=theta(2:end); %first and last point are the same. remove first point
if ~Edges %Edges==False
    psi=linspace(0,pi-Half*pi/2,Npsi+2);
    psi=psi(2:end-1); %remove edges. produced points will be redundent (as many as Npsi*2)
else
    psi=linspace(0,pi-Half*pi/2,Npsi);
end

%construct mesh
[Theta,Psi]=meshgrid(theta,psi);

%compute vertices
X=R*cos(Theta).*sin(Psi);
Y=R*sin(Theta).*sin(Psi);
Z=R*cos(Psi);

if exist('Trans','var') && ~isempty(Trans)
    X=X+Trans(1); Y=Y+Trans(2); Z=Z+Trans(3);
end
xyz=cat(3,X,Y,Z);
end
function [xyz,Z,Theta]=CalcCylinder(R,L,Nz,Ntheta,varargin)
%{
%INPUT:
%R - raidus
%Ntheta - amount of theta points (descrization). theta ranges [0,2pi]
%Nz -  amount of z points (descrization). z ranges [0,L]

%Varagin:
%Trans - translation [x,y,z]

%OUTPUT:
%xyz - size Nz x Ntheta x 3 where the depth dimension is constructed as [x,y,z]
%Theta, Z. matrices of size Ntheta x Nz
%}

%Obtain inputs
for ind=1:2:length(varargin)
    comm=lower(varargin{ind});
    switch comm
        case 'trans'
            Trans=varargin{ind+1};
    end
end

%create parameter vectors
theta=linspace(0,2*pi,Ntheta+1);
theta=theta(2:end); %first and last point are the same. remove first point
z=linspace(0,L,Nz);

%construct mesh
[Theta,Z]=meshgrid(theta,z);

%compute vertices
X=R*cos(Theta);
Y=R*sin(Theta);

if exist('Trans','var') && ~isempty(Trans)
    X=X+Trans(1); Y=Y+Trans(2); Z=Z+Trans(3);
end
xyz=cat(3,X,Y,Z);
end
function r=FindSphereRadius(xyz)
%{
%Assumption:
%xyz "spherical" half dome point cloud is centered around the origin facing
%up

%Algorithm:
%project point cloud on ZY and ZX plane
%calculate convex hall of both projections
%solve algebric sphereical fit (non linear)
%r=mean([rzx,rzy])

%Input
%xyz - matrix of [x,y,z] point cloud ((mx3)

%Output:
%r - radius of fitted sphere
%}
zy=xyz(:,[2,3]);
zyhull=zy(convhull(zy),:);
[~,rzy]=fitcircle(zyhull,'linear');

zx=xyz(:,[1,3]);
zxhull=zy(convhull(zx),:);
[~,rzx]=fitcircle(zxhull,'linear');

r=mean([rzx,rzy]);
end
function BlownPntCld=RadialNonRigidRegistration(Method,MovPntCld,StaticPntCld,Rbounds,Xcenter,varargin)
%{
%Blow up MovPntCld to StaticPntCld using fmincon (optimization) on
%Hausdorff distance to register MovPntCld onto StaticPntCld.
%Algorithm works in two steps:
%--rough registration to fit all MovPntCld with the same radius of expansion
%--a fine registration with fits different radii to different points.

%two methods for radial registration in this code:
%--cylinderical: radial registration from the z axis
%--spherical: radial registration from a given point

%Algorithm Assumption:
%--MovPntCld is bounded within StaticPntCld.
%--to register MovPntCld onto StaticPntCLd one must expand MovPntCld radially
%from a given point/axis.

%Example for expansion from a given axis:
%|  <-|   |->  |
%|  <-|   |->  |
%|  <-|   |->  |

%Input:
%Method - 'Cylinderical'/'Spherical'
%StaticPntCld - Point cloud in format of [x,y,z] (mx3) OR pointCloud
%object with similar .Location value
%MovPntCld - Point cloud in format mx3 or mxnx3 where the 3~[x,y,z]
%OR pointCloud object with similar.Location value
%RBounds - [Minimal radius,Maximal radius] to add to a point in MovPntCld
%Xcenter - [x,y,z] mx3 point in 3D space. if spherical mehthod was chosen it
%is the point from which radial expansion happens about. if cylinderical
%method was chosen, Xcenter(1,2) is the [x,y] coordante where the centeral axis passes
%through and expansion happens about

%Varargin inputs:
%DiffMinChange/DiffMaxChange - double, scalar. min/max change in radial
%values of points
%MaxIterations - integer,scalar.
%Display - 'none'/'iter/'final'. displays data in fmincon. see fmincon
%options for more
%OutputClass - 'double'/'pointCloud'. default set to 'double'

%Output:
%BlownPntCld - point cloud size format as MovPntCld(input). class may vary
%depending on Varargin input Outputclass choice.
%}

%varargin input and their defaults
DeltaR=Rbounds(2)-Rbounds(1);
DiffMinChange=DeltaR/1000; DiffMaxChange=DeltaR/10; MaxIterations=20;
Display='none'; OutputClass='double';
for ind=1:2:length(varargin)
    comm=lower(varargin{ind});
    switch comm
        case 'diffminchange'
            DiffMinChange=varargin{ind+1};
        case 'diffmaxchange'
            DiffMaxChange=varargin{ind+1};
        case 'maxiterations'
            MaxIterations=varargin{ind+1};
        case 'display'
            Display=varargin{ind+1};
        case 'outputclass'
            OutputClass=varargin{ind+1};
    end
end
%Default variables
if strcmpi(Method,'Spherical') && (~exist('Xcenter','var') || isempty(Xcenter))
    error('You must include Xcenter if spherical method was chosen. See varargin inputs');
end

%if point clouds give in uncomfi formats, change them
%StaticPntCld needs to be pointCloud object for cost function (use of
%findNearestNeighbors)
%MovPntCld should be in format mx3 or mxnx3 where the 3~[x,y,z]
if isa(MovPntCld,'pointCloud'), MovPntCld=MovPntCld.Location; end
if ~isa(StaticPntCld,'pointCloud'), StaticPntCld=pointCloud(StaticPntCld); end

%reshape MovPntCld to mx3 and save original size to reshape back to it if
%required
sz=size(MovPntCld);
if numel(sz)>2, MovPntCld=reshape(MovPntCld,[],3); end
m=size(MovPntCld,1); %find amount of points to move

%find t - matrix mx3 of radial direction per point
switch lower(Method)
    case 'spherical'
        t=(MovPntCld-Xcenter)./vecnorm(MovPntCld-Xcenter,2,2);
    case 'cylinderical'
        xy=MovPntCld(:,[1,2]);
        t=[xy-Xcenter(1:2),zeros(m,1)]./vecnorm(xy-Xcenter(1:2),2,2);
end

fun=@(r) CostFcn(MovPntCld,StaticPntCld,t,r);
options = optimoptions('fmincon','Algorithm','interior-point','Display',Display,...
    'DiffMaxChange',DiffMaxChange,'DiffMinChange',DiffMinChange,...
    'MaxIterations',MaxIterations,'MaxFunctionEvaluations',MaxIterations*m);
%Rough registration
lb=Rbounds(1);
ub=Rbounds(2);
r0=DiffMinChange;
rRough=fmincon(fun,r0,[],[],[],[],lb,ub,[],options);
%Fine registration
lb=Rbounds(1)*ones(m,1);
ub=Rbounds(2)*ones(m,1);
r0=rRough*ones(m,1);
rFine=fmincon(fun,r0,[],[],[],[],lb,ub,[],options);

%Create output
BlownPntCld=MovPntCld+t.*rFine;
if numel(sz)>2, BlownPntCld=reshape(BlownPntCld,sz); end
if strcmpi(OutputClass,'pointCloud'), BlownPntCld=pointCloud(BlownPntCld); end

    function Cost=CostFcn(MovPntCld,StaticPntCld,t,r)
        %note: Nested functions variables who share a common name with
        %their parent function will become global. for such reason, 'N' was
        %used instead of 'm'
        N=size(MovPntCld,1);
        Cost=0;
        points=MovPntCld+t.*r;
        for i=1:N
            [~,d]=findNearestNeighbors(StaticPntCld,points(i,:),1);
            Cost=Cost+d^2;
        end
    end
end