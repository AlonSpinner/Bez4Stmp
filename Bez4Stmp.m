classdef Bez4Stmp
    %For final project in faculty of mechanichal engineering, Technion
    %2019.
    %By Yam Ben Natan and Alon Spinner
    properties (SetAccess = private)  %equivalent to "read only" (Getable but not Setable)
        Scan %computed in constructor. Point Cloud provided to algorithm
        PointCloud %computed in constructor. PointCloud after processing (Transform2Z and downsampling)
        BoundingCylinderLength %computed in constructor
        BoundingCylinderRadius %computed in constructor
        Xnode %computed in constructor
        CompactControlPoints %control points in initial form (before RadialNonRigid registration)
        BlownControlPoints %control points after Compact has been blown up
    end
    properties (SetAccess = public) %can be changed by user
        GridLengthFcn=@(R,L) R*(0.5*R/L); %after downsampling: VerticesAmount*(GridLength^2)/(2*pi*R*L)~1 by construction
        PatchAmount=20; %default
        BezierOrder=3; %default
        XnodeCalculationMethod='normalSTD'; %default
    end
    methods %public methods, accessible from outside the class
        function obj=Bez4Stmp(PntCld) %Constructor
            %Input:
            %PntCld - points to xyz double matrix of size mx3 point
            %cloud OR pointCloud class object.
            
            %depending on it's class, assumes the following:
            %double: xyz double matrix of size mx3 [x,y,z]
            %char: file path
            
            %Output:
            %initalized class
            
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
            switch class(PntCld)
                case 'double'
                    xyz=PntCld;
                case 'char'
                    FileName=PntCld;
                    [~,~,ext]=fileparts(FileName);
                    switch ext
                        case '.stl'
                            [~,xyz]=stlread(FileName); %obtain vertices of STL
                        case '.mat'
                            load(FileName,'xyz'); %assumes mat file contains variable named xyz and loads it
                            if isa(xyz,'pointCloud'), xyz=xyz.Location; end %if point cloud and not mx3 double. make it mx3 double
                    end
                case 'pointCloud' %incase user went ahead and created a point cloud object
                    xyz=PntCld.Location;
                otherwise
                    clear obj %delete default object which was created when method was called
                    error(sprintf(['Wrong Input\n\n',...
                        'Input must be one of the following:\n',...
                        ' -mx3 double matrix containing [x,y,z]\n',...
                        ' -pointCloud class variable\n',...
                        ' -file name of stl\n',...
                        ' -file name of .m file containing pointCloud or mx3 double matrix named "xyz"']));
            end
            %% Run
            %introduce to parameters before any preprocessing to object
            obj.Scan=pointCloud(xyz);
            %align with z axis, sets base on z=0, sort vertices by z value
            T2Z_xyz=Transform2Z(xyz);
            %find Radius and Length of bounding Cylinder
            [R,L]=BoundingCylinder(T2Z_xyz);
            obj.BoundingCylinderLength=L; %introduce to parameters
            obj.BoundingCylinderRadius=R; %introduce to parameters
            
            %Create PointCloud - post preprocessing point cloud object
            PtCld=pointCloud(T2Z_xyz); %create point cloud object from preprocessed point cloud
            GridLength=obj.GridLengthFcn(R,L); %obtain voxel edge length by default function
            PtCld=pcdownsample(PtCld,'gridAverage',GridLength); %Downsample by grid legnth
            obj.PointCloud=PtCld; %introduce to parameters
            
            %Find Xnode - point seperating between sphereical and
            %Cylinderical part
            Xnode=FindXnode(PtCld,L,obj.XnodeCalculationMethod);
            obj.Xnode=Xnode; %introduce to parameters
            
            %Create Compact control points - control points before blowing
            %algorithm
            P=BezBellCurve(R/2,R/2,Xnode(3)+R/2,20);
            xyz=RotateXZcurve(P,20);
            xyz(:,:,1)=xyz(:,:,1)+Xnode(1); xyz(:,:,2)=xyz(:,:,2)+Xnode(2);
            obj.CompactControlPoints=pointCloud(xyz);
            
            %Create BlownControlPoints - control points after blowing
            %algorithm
            h=helpdlg(sprintf(['Creating BlownControlPoints from CompactControlPoints.',...
                'This might take a few mintues']));
            %Split process to two parts - sphere and cylinder. then merge
            %Sphere:
            roiSph=[-R,R,-R,R,Xnode(3),L]; %range of intereset (bounding box [xbound,ybound,zbound])
            StatSph=select(obj.PointCloud,findPointsInROI(obj.PointCloud,roiSph)); %stationary sphere from obj.PointCloud
            MovSph=select(obj.CompactControlPoints,findPointsInROI(obj.CompactControlPoints,roiSph)); %sphere to blow (from obj.CompactControlPoints)
            BlownSph=RadialNonRigidRegistration('Spherical',MovSph,StatSph,L-Xnode(3),'Xnode',Xnode,'Display','iter'); %blown sphere to be  merged
            %Cylinder:
            roiCyl=[-R,R,-R,R,0,Xnode(3)]; %range of intereset (bounding box [xbound,ybound,zbound])
            StatCyl=select(obj.PointCloud,findPointsInROI(obj.PointCloud,roiCyl)); %stationary cylinder from obj.PointCloud
            MovCyl=select(obj.CompactControlPoints,findPointsInROI(obj.CompactControlPoints,roiCyl)); %cylinder to blow (from obj.CompactControlPoints)
            BlownCyl=RadialNonRigidRegistration('Cylinderical',MovCyl,StatCyl,R,'Display','iter'); %blown cylinder to be  merged
            %Merge:
            obj.BlownControlPoints=pcmerge(BlownSph,BlownCyl,GridLength); %merge
            close(h); %close helpdlg
        end
    end
    methods (Static)
        function Ax=DrawPointCloud(PntCld,varargin)
            %INPUT:
            %PointCloud - can be of differnent types:
            %numeric matrix of size mx3
            %numeric matrix of size mxnx3 where the third element
            %corresponds to [x,y,z]
            %pointCloud class object
            %varargin:
            %Color - color as [R,G,B] or matlab strings
            %ColorMap - colormap as string ('jet','parula'...)
            %Msize - numeric value for markersize
            %Ax - handle of axes
            %title - title of axes
            
            %OUTPUT:
            %CldPHandle - a handle to the axes object
            
            %check input. if is double, create a pointCloud object
            switch class(PntCld)
                case 'double'
                    sz=size(PntCld);
                    if ~(numel(sz)==2 && sz(2)==3) && ~(numel(sz)==3 && sz(3)==3)
                        error(sprintf(['When introducing "double" PntCld make sure',...
                            'its of size nx3 correorsponding to [x,y,z] or mxnx3',...
                            'where the third size element corresponds to [x,y,z]']));
                    end
                case 'pointCloud'
                    PntCld=PntCld.Location;
            end
            
            %Obtain varargin inputs
            for ind=1:2:length(varargin)
                comm=lower(varargin{ind});
                switch comm
                    case 'ax'
                        Ax=varargin{ind+1};
                    case 'colormap'
                        ColorMap=varargin{ind+1};
                    case 'color'
                        Color=varargin{ind+1};
                    case 'msize'
                        Msize=varargin{ind+1};
                    case 'title'
                        Title=varargin{ind+1};
                end
            end
            
            %varargin input check
            %if both Color and ColorMap are provided, only color will be relevant.
            if exist('ColorMap','var') && exist('Color','var')
                warning('both ColorMap and Color were provided. Draws by ColorMap');
            end
            
            %Default values and some more input check (in Color)
            if ~exist('Ax','var') || isempty(Ax)
                Fig=figure('color',[0,0,0]);
                Ax=axes(Fig,'color',[0,0,0],'ZColor',[1,1,1],'XColor',[1,1,1],'YColor',[1,1,1]);
                xlabel(Ax,'x'); ylabel(Ax,'y'); zlabel(Ax,'z');
                axis(Ax,'equal'); grid(Ax,'on'); hold(Ax,'on'); view(Ax,3);
            end
            if exist('Title','var') && isa(Title,'char') %set title if provided
                title(Ax,['\color{white}',Title]);
            end
            if ~exist('Msize','var') || isempty(Msize)
                Msize=1;
            end
            if ~exist('ColorMap','var') || isempty(ColorMap)
                ColorMap='parula';
            end
            
            %draw. function call pcshow depends on color choice
            if ~exist('Color','var') || isempty(Color)
                pcshow(PntCld,'MarkerSize',Msize,'Parent',Ax);
                colormap(Ax,ColorMap);
            else
                if numel(Color)==3
                    %         Color=Color.*ones(size(PntCld));
                    pcshow(PntCld,Color,'MarkerSize',Msize,'Parent',Ax);
                else
                    error('Color must be a 3x1 or 1x3 RGB vector');
                end
            end
        end
    end %methods that have no interaction with the object created by the class
end
%% algorithm pipeline functions
function Newxyz=Transform2Z(xyz)
%Translate xyz points so that CGxy=[0,0], floor  is z=0, and primary
%axis of the convex shape will be the Z axis.
%NewXYZ is also sorted by Z (first row with minimal Z and last row with
%maximal)

%Input:
%xyz - [x,y,z] numeric matrix mX3

%Output:
%Newxyz - [x,y,z] numeric matrix mX3.

%centerlize xyz and compute pca analysis
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
%Assumption:
%xyz was already processed by function Transfom2Z:
%it is aligned with z axis, has base on z=0

%Input
%xyz - matrix of [x,y,z] point cloud

%Output:
%R - radius of bounding cylinder
%L - Length of bouding cylinder

% Calculate L
L=max(xyz(:,3));

% Calculate R
% algerbric error minimization (not geometric fit)
%
% Reference: "Least-squares fitting of circles and ellipses", W. Gander,
% G. Golub, R. Strebel - BIT Numerical Mathematics, 1994, Springer

% This implementation copyright Richard Brown, 2007, but is freely
% available to copy, use, or modify as long as this line is maintained

% Obtained from "fitcircle.m" from mathworks file exchange
% https://www.mathworks.com/matlabcentral/fileexchange/15060-fitcircle-m

%project all XYZ points to XY plane and find convex hull.
%Circl will be  fitted on convex hull
xy=xyz(:,[1,2]);
CHull=xy(convhull(xy),:);  %[x,y] mat. find convex hull points

%Convenience variables
m =size(CHull,1);
x1=CHull(:,1);
x2=CHull(:,2);

% Compute best fit w.r.t. algebraic error using linear least squares
%
% Circle is represented as a matrix quadratic form
%     ax'x + b'x + c = 0
% Linear least squares estimate found by minimising Bu = 0 s.t. norm(u) = 1
%     where u = [a; b; c]

% Form the coefficient matrix
B = [x1.^2 + x2.^2, x1, x2, ones(m, 1)];

% Least squares estimate is right singular vector corresp. to smallest
% singular value of B
[~,~, V] = svd(B);
u = V(:, 4);

% For clarity, set the quadratic form variables
a = u(1);
b = u(2:3);
c = u(4);

% Convert to centre/radius
z = -b / (2*a);
R = sqrt((norm(b)/(2*a))^2 - c/a);
end
function Xnode=FindXnode(PtCld,L,method,varargin)
%Assumption:
%--PtCld was already processed by function Transfom2Z:
%it is aligned with z axis, has base on z=0
%--PtCld has already gone through downsampling using GridLengthFcn(R,L) and
%function pcdownsample.

%Input:
%PtCldP  - pointCloud object
%L - BoundingCylinderLength
%Method - 'GeomtreyAssumption'/'normalSTD'/'normalThreshold'
%Varargin - depending on method:
%normalSTD - BoundingCylinderLength
%GeomtreyAssumption: BoundingCylinderRadius
%normalThreshold: requires Threshold scalar ranges from [0,1]
%suggestion: 0.4

%Output:
%Xnode [x,y,z] is the node that seprates between the "sphere" part and the
%"cylinder" part in the given cloud point

%Assumptions:
%CldP main axis
%     ___
%   /     \
% /        |
%|     x   |
%|         |
%|         |
%|         |


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
        pass=find(PtCld.Location(:,3)>L/4); %all vertices that pass the condition
        CleanPtCld=select(PtCld,pass);
        k=8;
        normals=pcnormals(CleanPtCld,k);
        Idx=findchangepts(normals(:,3),'Statistic','std'); %check for the point in which the std changes the moves
        Znode=CleanPtCld.Location(Idx,3);
    case 'normalThreshold'
        if ~exist('Threshold','var') || isempty(Threshold) || (Threshold>1 || Threshold<0)
            error('Failed to compute Xnode. Need to provide a Threshold in range [0,1] with the chosen method')
        end
        %sometimes there is noise in the lower parts of the given scan that prodces vertical normals
        pass=find(PtCld.Location(:,3)>L/4); %all vertices that pass the condition
        CleanPtCld=select(PtCld,pass);
        k=8; %see explination in normalSTD method.
        normals=pcnormals(CleanPtCld,k);
        pass=find(abs(normals)>Threshold); %all vertices that pass the condition
        Idx=pass(1); %take first vertex that pass condition
        Znode=CleanPtCld.Location(Idx,3);
end

%find [x,y] of Xnode by mean on all the points above it (higher than Znode)
pass=find(PtCld.Location(:,3)>Znode); %all vertices that pass the condition
AbvPtCld=select(PtCld,pass);
Xnode=[mean(AbvPtCld.Location(:,[1,2])),Znode];

end
function BlownPntCld=RadialNonRigidRegistration(Method,MovPntCld,StaticPntCld,MaxRadius,varargin)
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
%MovPntCld - Point cloud cloud in format mx3 or mxnx3 where the 3~[x,y,z]
%OR pointCloud object with similar.Location value
%MaxRadius - Maximal radius to add to a point in MovPntCld
%Xnode - point from with

%Varargin inputs:
%Xnode - [x,y,z] mx3 point in 3D space. if spherical mehthod was chosen it
%is the point from which radial expansion happens about
%DiffMinChange/DiffMaxChange - double, scalar. min/max change in radial
%values of points
%MaxIterations - integer,scalar.
%Display - 'none'/'iter/'final'. displays data in fmincon. see fmincon
%options for more

%Output:
%BlownPntCld - pointCloud object with .Location in the same size format as
%MovPntCld(input)

%varargin input and their defaults
DiffMinChange=MaxRadius/1000; DiffMaxChange=MaxRadius/10; MaxIterations=20;
Display='none';
for ind=1:2:length(varargin)
    comm=lower(varargin{ind});
    switch comm
        case 'xnode'
            Xnode=varargin{ind+1};
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
%Default variables
if strcmpi(Method,'Spherical') && (~exist('Xnode','var') || isempty(Xnode))
    error('You must include Xnode if spherical method was chosen. See varargin inputs');
end


%Obtain double matrices if pointCloud objects were given
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
        t=(MovPntCld-Xnode)./vecnorm(MovPntCld-Xnode,2,2);
    case 'cylinderical'
        xy=MovPntCld(:,[1,2]);
        t=[xy,zeros(m,1)]./vecnorm(xy,2,2);
end

fun=@(r) CostFcn(MovPntCld,StaticPntCld,t,r);
options = optimoptions('fmincon','Algorithm','interior-point','Display',Display,...
    'DiffMaxChange',DiffMaxChange,'DiffMinChange',DiffMinChange,...
    'MaxIterations',MaxIterations);
%Rough registration
lb=0;
ub=MaxRadius;
r0=DiffMinChange;
rRough=fmincon(fun,r0,[],[],[],[],lb,ub,[],options);
%Fine registration
lb=zeros(m,1);
ub=MaxRadius*ones(m,1);
r0=rRough*ones(m,1);
rFine=fmincon(fun,r0,[],[],[],[],lb,ub,[],options);

%Create output
MovPntCld=MovPntCld+t.*rFine;
if numel(sz)>2, MovPntCld=reshape(MovPntCld,sz); end
BlownPntCld=pointCloud(MovPntCld);

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
%% extra functions (note: can't call one method from another)
function varargout = stlread(file)
%Works only with binary STLS
%imported from https://www.mathworks.com/matlabcentral/fileexchange/22409-stl-file-reader
%slightly edited to remove ASCII parts which did not work
%anyhow

% STLREAD imports geometry from an STL file into MATLAB.
%    FV = STLREAD(FILENAME) imports triangular faces from the ASCII or binary
%    STL file idicated by FILENAME, and returns the patch struct FV, with fields
%    'faces' and 'vertices'.
%
%    [F,V] = STLREAD(FILENAME) returns the faces F and vertices V separately.
%
%    [F,V,N] = STLREAD(FILENAME) also returns the face normal vectors.
%
%    The faces and vertices are arranged in the format used by the PATCH plot
%    object.
% Copyright 2011 The MathWorks, Inc.
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
function [F,V,N] = stlbinary(M)
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
function [P,Q]=BezBellCurve(R,d,L,Nq)
%Builds a 2D bezier curve and evalutes Nq points on it
%(eqidistance on parameter q in [0,1])

%    Control points (symbolised by "@") are set as such:
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
%Q - Control points [x,y,z] (mx2)

%Example:
% L=10; R=3; d=3; Nq=10;
% [P,Q]=Bez4Stmp.BellCurve(R,d,L,Nq);
% x=axes; hold(x,'on'); grid(x,'on'); axis(x,'equal');
% scatter(P(:,1),P(:,2),5,'b');
% scatter(Q(:,1),Q(:,2),10,'r');
% %NOTE: d=R, d=0 and d=L are the three interesting cases

%build control points and parameter vector
Q=[R,0;
    R,L-d;
    R,L;
    0,L];
q=linspace(0,1,Nq);

%Evaluate points
P=zeros(Nq,2);
for i=1:Nq
    P(i,:)=BezCrv_DeCasteljau(Q,q(i));
end

end
function xyz=RotateXZcurve(xz,n)
%Input:
%xz [x,z] of evalualted points on 2D curve (mx2 array)
%n - multiplyer amount

%Output:
%xyz  evaluated points on 3D surface (mxn3 array) where depth dimension is
%       orginized [x,y]z.
%       created by turnning the xy curve n times around the y axis. each
%       time by 2*pi/N amount

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
function R=BezCrv_DeCasteljau(Q,q)
%Evaluates Bezier Curve by given control points and parameter
%value

%Q - control points in format [x] or [x,y,z] dimension matrix. top row is q=0.
%q - running parameter of bezier curve. 0=<q<=1
%R - [x] or [x,y,z] format. point on bezier curve
% https://pages.mtu.edu/~shene/COURSES/cs3621/NOTES/spline/Bezier/de-casteljau.html

n=size(Q,1)-1; %degree of bezier polynomial
for k=1:(n+1)
    for i=1:(n+1-k) %doesnt enter in first iteration
        Q(i,:)=(1-q)*Q(i,:)+q*Q(i+1,:);
    end
end
R=Q(1,:);
end
%% functions not in use
function [xyz,Theta,Psi]=CalcSphere(R,Ntheta,Npsi,varargin)
%INPUT:
%R - raidus
%Ntheta - amount of theta points (descrization). theta ranges [0,2pi]
%Npsi - amount of psi points (descrization). psi ranges [0,pi]

%Varagin:
%Trans - translation [x,y,z]
%Half - true/false. true - plots half dome (pointing up)

%OUTPUT:
%xyz - size Ntheta x Npsi x 3 where the depth dimension is constructed as [x,y,z]
%Theta, Psi. matrices of size Ntheta x Npsi

%default values:
Half=false;
%Obtain inputs
for ind=1:2:length(varargin)
    comm=lower(varargin{ind});
    switch comm
        case 'trans'
            Trans=varargin{ind+1};
        case 'half'
            Half=varargin{ind+1};
    end
end

%create parameter vectors
theta=linspace(0,2*pi,Ntheta+1);
theta=theta(2:end); %first and last point are the same. remove first point
psi=linspace(0,pi-Half*pi/2,Npsi+2);
psi=psi(2:end-1); %remove edges. produced points will be redundent (as many as Npsi*2)

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
function [xyz,Theta,Z]=CalcCylinder(R,L,Ntheta,Nz,varargin)
%INPUT:
%R - raidus
%Ntheta - amount of theta points (descrization). theta ranges [0,2pi]
%Nz -  amount of z points (descrization). z ranges [0,L]

%Varagin:
%Trans - translation [x,y,z]

%OUTPUT:
%xyz - size Ntheta x Nz x 3 where the depth dimension is constructed as [x,y,z]
%Theta, Z. matrices of size Ntheta x Nz

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
function P=CircleBellCrv(R,h,Ncyl,Nsph)
%Input:
%R - radius
%h - height
%Ncyl - number of points in "cylinderical" part
%Nsph - number of points in "Spherical" part

%Output:
%P - [x,y] (size (Nsph+Ncyl)x2 of evaluated points

%---           ^
%  ^  \        |
%  |    \      |
%  R      \    |
%  |       \   |
%  v       |   h
%          |   |
%          |   |
% <---R--->|   v

%built function. note: ((z-(h-R))/R)=sin(theta)
fcnrz=@(z) R-heaviside(z-(h-R))*R*(1-sqrt(1-((z-(h-R))/R)^2));

%build vectors
zcyl=linspace(0,h-R/2,Ncyl);
zsph=linspace(h-R/2,h,Nsph+1);
zsph=zsph(2:end); %first point was as zcyl last point
z=[zcyl,zsph];

%evalulate points
r=zeros(size(z));
for i=1:N
    r(i)=fcnrz(z(i));
end

P=[r(:),z(:)];
end