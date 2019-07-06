classdef Bez4Stmp
    %For final project in faculty of mechanichal engineering, Technion
    %2019.
    %By Yam Ben Natan and Alon Spinner
    properties (SetAccess = private)  %equivalent to "read only" (Getable but not Setable)
        %all properties below are initalized in constructor
        
        Scan %pointCloud mx3. Point Cloud provided to algorithm
        PointCloud %pointCloud mx3. PointCloud after processing (Transform2Z and downsampling)
        BoundingCylinderLength %%double 1x1.
        BoundingCylinderRadius %%double 1x1.
        SphereRadius %double 1x1.
        Xcenter %double [x,y,z].
        CompactNodes %pointCloud mxnx3. Nodes in initial form (before RadialNonRigid registration)
        BlownNodes %pointCloud mxnx3. nodes after CompactNodes have been thorugh RadialNonRigid registration
        CPcyl %struct of bezier mesh control points containing info of vertices, patches and connectivity
        CPsph
    end
    properties (SetAccess = public) %can be changed by user
        GridLengthFcn=@(R,L) R*(0.5*R/L); %after downsampling: VerticesAmount*(GridLength^2)/(2*pi*R*L)~1 by construction
        BezierPatchAmount=16; %default
        SphLayers=2; %default
        CylLayers=2;
        Slices=4;
        BezierOrder=3; %default
        XcenterCalculationMethod='normalSTD'; %default
    end
    methods %public methods, accessible from outside the class
        function obj=Bez4Stmp(PntCldIn) %Constructor
            %Input:
            %PntCldIn - points to xyz double matrix of size mx3 point
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
            switch class(PntCldIn)
                case 'double'
                    xyz=PntCldIn;
                case 'char'
                    FileName=PntCldIn;
                    [~,~,ext]=fileparts(FileName);
                    switch ext
                        case '.stl'
                            [~,xyz]=stlread(FileName); %obtain vertices of STL
                        case '.mat'
                            load(FileName,'xyz'); %assumes mat file contains variable named xyz and loads it
                            if isa(xyz,'pointCloud'), xyz=xyz.Location; end %if point cloud and not mx3 double. make it mx3 double
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
            PntCld=pointCloud(T2Z_xyz); %create point cloud object from preprocessed point cloud
            GridLength=obj.GridLengthFcn(R,L); %obtain voxel edge length by default function
            PntCld=pcdownsample(PntCld,'gridAverage',GridLength); %Downsample by grid legnth
            obj.PointCloud=PntCld; %introduce to parameters
            
            %Find Xcenter - point seperating between sphereical and
            %Cylinderical part
            Xcntr=FindXcenter(PntCld,L,obj.XcenterCalculationMethod);
            obj.Xcenter=Xcntr; %introduce to parameters
            
            %Find radius of sphereical part
            r=L-Xcntr(:,3);
            obj.SphereRadius=r;
            
            %Create Compact nodes - nodes before blowing
            %algorithm
            [Ncyl,Nsph,Ncircum]=MeshNodesAmount(obj.SphLayers,obj.CylLayers,obj.Slices,obj.BezierOrder);
            SphNodes=CalcSphere(R,Nsph,Ncircum,'Trans',Xcntr,'half',true);
            CPsph=StumpMesh2CP(SphNodes,obj.BezierOrder);
            CylNodes=CalcCylinder(R,Xcntr(3),Ncyl,Ncircum);
            CPcyl=StumpMesh2CP(CylNodes,obj.BezierOrder);
            
%             CompactSphInd=(P(:,2)>Xcntr(3)); CompactCylInd=~CompactSphInd; %boolean arrays [0,0,0,1,1,1,1]
%             xyz=RotateXZcurve(P,Ncircum); %rotate to create double size Nvert x Ncircum x 3
%             xyz(:,:,1)=xyz(:,:,1)+Xcntr(1); xyz(:,:,2)=xyz(:,:,2)+Xcntr(2); %align vertical axis to pass through Xcenter
%             obj.CompactNodes=pointCloud(xyz); %introduce to object parameters
%             CompactSph=xyz(CompactSphInd,:,:); CompactCyl=xyz(CompactCylInd,:,:); %seperate for expansion(blowing) algorithm
            
            %Create BlownNodes - nodes after blowing algorithm
            h=helpdlg(sprintf(['Creating BlownNodes from CompactNodes.\n',...
                'This might take a few mintues']));
            %Split process to two parts - sphere and cylinder. then merge
            StatSphInd=PntCld.Location(:,3)>Xcntr(3); StatCylInd=~StatSphInd;%boolean arrays [0,0,0,1,1,1,1]
            %Sphere:
            StatSph=PntCld.Location(StatSphInd,:);
            CPsph=RadialOptimization4CP('Spherical',CPsph,StatSph,[-r/2,2*r],Xcntr,'Display','iter'); %blown sphere to be  merged
            %Cylinder:
            StatCyl=PntCld.Location(StatCylInd,:);
            CPcyl=RadialOptimization4CP('Cylinderical',CPcyl,StatCyl,[-R/2,2*R],Xcntr,'Display','iter'); 
            %Merge:
            close(h); %close helpdlg
            
            %create CP struct and introduce to parameters
            obj.CPcyl=CPcyl;
            obj.CPsph=CPsph;
        end
        function Handles=DrawBezierPatches(obj,varargin)
            %draws bezier parameteric surfaces from data stored in obj.CP
            
            %Varagin inputs:
            %Ax - handle of axes. if none chosen will draw to gca
            %N - number of points drawn are N^2. by default N=30
            %Curvature - 'gaussian'/'mean'/'none'. default set to 'none'
            %Color - color of all patches. if not specified - will be by lines (matlab colormap).
            %PauseTime - pause time between patch drawings (seconds).
                    %default set to "none"
              
            %note: Providing both Color and Curvature returns an error.
            
            %Outputs:
            %Handles - array of handles size(1,Patch Amount)
            
            %Note:
            %CP is a struct with:
            
            %CP.v - vertices [X,Y,Z] matrix format
            
            %CP.p -  patches. size([BezO+1,BezO+1,Patch Amount]).
            %Values are row index of points in CP.v. depth index is patch indexing.
            
            %CP.c - patch connectivity. size([Ptch Amount,4]). The row index of CP.c indicates the patch involved.
            %Values are patch indexes. a patch can have 4 neighbors and they are ordered [Left, Right, Up,
            %Down] in the row. if a patch has no connectivity in a certin direction,
            %the related value will be 0.
            
            pCP=obj.CPcyl; %pCP - private CP (just to shorten writing without referencing from obj)
            
            %default values
            N=30; Curvature='none'; PauseTime=0;
            %Obtain inputs
            for ind=1:2:length(varargin)
                comm=lower(varargin{ind});
                switch comm
                    case 'n'
                        N=varargin{ind+1};
                    case 'ax'
                        Ax=varargin{ind+1};
                    case 'color'
                        Color=varargin{ind+1};
                    case 'curvature'
                        Curvature=varargin{ind+1};
                    case 'pausetime'
                        PauseTime=varargin{ind+1};
                end
            end
            
            %Check input validity and create color matrix if required (==if
            %drawing is NOT by curvature)
            PtchAmnt=size(pCP.p,3);
            if strcmpi(Curvature,'none')  %if Curvature wasnt provided or set to "none"
                if ~exist('Color','var') || isempty(Color) %if color wasnt provided
                    PtchColor=lines(PtchAmnt);
                else
                    PtchColor=repmat(Color,PtchAmnt,1); %if color was provided
                end
            else %if Curvature was provided and set to "gaussian"/"mean"
                if exist('Color','var') %if color exists, return error
                    error('Color can only be provided if curvature is set to "none"');
                end
            end
            
            %create axes if not given
            if ~exist('Ax','var') || isempty(Ax)
                Fig=figure('color',[0,0,0]);
                Ax=axes(Fig,'color',[0,0,0],'ZColor',[1,1,1],'XColor',[1,1,1],'YColor',[1,1,1]);
                xlabel(Ax,'x'); ylabel(Ax,'y'); zlabel(Ax,'z');
                axis(Ax,'equal'); grid(Ax,'on'); hold(Ax,'on'); view(Ax,3);
            end
            
            %initalization for drawing patches 
            Handles=gobjects(1,PtchAmnt); %initalize handle array
            uOp1=size(pCP.p,2); %u order plus 1 (amount of control points in the u direction)
            vOp1=size(pCP.p,1);
            
            %Draw patches       
            for k=1:PtchAmnt
                %Evaluate points of patch
                Pcp=reshape(pCP.v(pCP.p(:,:,k),:),[vOp1,uOp1,3]); %Obtain patch control points. format explained above
                Psrfp=BezPtchCP2SrfP(Pcp,N); %compute surface points
                x=Psrfp(:,:,1); y=Psrfp(:,:,2); z=Psrfp(:,:,3);
               
                %Draw evaluated points with surf command
                %Note: switch coded in forloop to allow cool "build up" graphics
                        %when plotting
                switch Curvature
                    case 'none'
                        Handles(k)=surf(Ax,x,y,z,...
                            'facecolor',PtchColor(k,:),'edgecolor','none','facealpha',0.6,...
                            'UserData',k);
                    case 'gaussian'
                        K=SurfCurvature(x,y,z);
                        K=filloutliers(K,'nearest');
                        Handles(k)=surf(Ax,x,y,z,K,'facecolor','interp','edgecolor','none');
                    case 'mean'
                        [~,H]=SurfCurvature(x,y,z);
                        Handles(k)=surf(Ax,x,y,z,H,'facecolor','interp','edgecolor','none');
                end               
                pause(PauseTime);
            end
            
            
            
            pCP=obj.CPsph; %pCP - private CP (just to shorten writing without referencing from obj)
            
            %default values
            N=30; Curvature='none'; PauseTime=0;
            %Obtain inputs
            for ind=1:2:length(varargin)
                comm=lower(varargin{ind});
                switch comm
                    case 'n'
                        N=varargin{ind+1};
                    case 'ax'
                        Ax=varargin{ind+1};
                    case 'color'
                        Color=varargin{ind+1};
                    case 'curvature'
                        Curvature=varargin{ind+1};
                    case 'pausetime'
                        PauseTime=varargin{ind+1};
                end
            end
            
            %Check input validity and create color matrix if required (==if
            %drawing is NOT by curvature)
            PtchAmnt=size(pCP.p,3);
            if strcmpi(Curvature,'none')  %if Curvature wasnt provided or set to "none"
                if ~exist('Color','var') || isempty(Color) %if color wasnt provided
                    PtchColor=lines(PtchAmnt);
                else
                    PtchColor=repmat(Color,PtchAmnt,1); %if color was provided
                end
            else %if Curvature was provided and set to "gaussian"/"mean"
                if exist('Color','var') %if color exists, return error
                    error('Color can only be provided if curvature is set to "none"');
                end
            end
            
            %create axes if not given
            if ~exist('Ax','var') || isempty(Ax)
                Fig=figure('color',[0,0,0]);
                Ax=axes(Fig,'color',[0,0,0],'ZColor',[1,1,1],'XColor',[1,1,1],'YColor',[1,1,1]);
                xlabel(Ax,'x'); ylabel(Ax,'y'); zlabel(Ax,'z');
                axis(Ax,'equal'); grid(Ax,'on'); hold(Ax,'on'); view(Ax,3);
            end
            
            %initalization for drawing patches 
            Handles=gobjects(1,PtchAmnt); %initalize handle array
            uOp1=size(pCP.p,2); %u order plus 1 (amount of control points in the u direction)
            vOp1=size(pCP.p,1);
            
            %Draw patches       
            for k=1:PtchAmnt
                %Evaluate points of patch
                Pcp=reshape(pCP.v(pCP.p(:,:,k),:),[vOp1,uOp1,3]); %Obtain patch control points. format explained above
                Psrfp=BezPtchCP2SrfP(Pcp,N); %compute surface points
                x=Psrfp(:,:,1); y=Psrfp(:,:,2); z=Psrfp(:,:,3);
               
                %Draw evaluated points with surf command
                %Note: switch coded in forloop to allow cool "build up" graphics
                        %when plotting
                switch Curvature
                    case 'none'
                        Handles(k)=surf(Ax,x,y,z,...
                            'facecolor',PtchColor(k,:),'edgecolor','none','facealpha',0.6,...
                            'UserData',k);
                    case 'gaussian'
                        K=SurfCurvature(x,y,z);
                        K=filloutliers(K,'nearest');
                        Handles(k)=surf(Ax,x,y,z,K,'facecolor','interp','edgecolor','none');
                    case 'mean'
                        [~,H]=SurfCurvature(x,y,z);
                        Handles(k)=surf(Ax,x,y,z,H,'facecolor','interp','edgecolor','none');
                end               
                pause(PauseTime);
            end
            
            
            
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
        function [hd,pInd,qInd]=Hausdorff(P,Q,varargin)
            % Calculates the Hausdorff Distance, hd, between two sets of points, P and
            % Q (which could be two trajectories). Sets P and Q must be matrices with
            % an equal number of columns (dimensions), though not necessarily an equal
            % number of rows (observations).
            
            %Input:
            %P - double matrix m1xn
            %Q - double matrix m2xn
            
            %varargin Input:
            %'OneSide' - loops on only on P for minimal distances
            
            %Output:
            %hd -scalar hausdorff distance
            %pInd,qInd - hd=norm(P(pIind,:)-Q(qInd,:))
            
            Method='TwoSide'; %classic hausdorff defintion
            if numel(varargin)>0 && strcmpi(varargin{1},'OneSide')
                Method='OneSide';
            end
            
            szP=size(P); szQ=size(Q);
            if szP(2)~=szQ(2), error('dimension size (number of columns) of both inputs need to be the same'); end
            
            ed2=@(p,q) sum((p-q).^2,2); %euclidian distant squared
            [pInd,qInd]=deal(ones(2,1));
            hd2P=0; %initalize
            
            %loop on P
            for ip=1:szP(1)
                [hd2,iq]=min(ed2(P(ip,:),Q)); %hausdorff distance
                if hd2>hd2P
                    hd2P=hd2;
                    pInd(1)=ip;
                    qInd(1)=iq;
                end
            end
            %loop on Q
            hd2Q=0;
            if strcmp(Method,'TwoSide')
                for iq=1:szQ(1)
                    [hd2,ip]=min(ed2(Q(iq,:),P)); %hausdorff distance
                    if hd2>hd2Q
                        hd2Q=hd2;
                        pInd(2)=ip;
                        qInd(2)=iq;
                    end
                end
            end
            %conclude
            [hd2,Ind]=max([hd2P,hd2Q]);
            pInd=pInd(Ind);
            qInd=qInd(Ind);
            hd=sqrt(hd2);
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
%xyz - matrix of [x,y,z] point cloud ((mx3)

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

%project all XYZ points to XY plane and find convex hull.
%Circl will be  fitted on convex hull
xy=xyz(:,[1,2]);
CHull=xy(convhull(xy),:);  %[x,y] mat. find convex hull points

[~,R]=fitcircle(CHull,'linear');
end
function Xcenter=FindXcenter(PtCld,L,method,varargin)
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
%Xcenter [x,y,z] is the node that seprates between the "sphere" part and the
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
            error('Failed to compute Xcenter. Need to provide a Threshold in range [0,1] with the chosen method')
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

%find [x,y] of Xcenter by mean on all the points above it (higher than Znode)
pass=find(PtCld.Location(:,3)>Znode); %all vertices that pass the condition
AbvPtCld=select(PtCld,pass);
Xcenter=[mean(AbvPtCld.Location(:,[1,2])),Znode];

end
function [Ncyl,Nsph,Ncircum]=MeshNodesAmount(SphLayers,CylLayers,Slices,BezierOrder)
%Input:
%CylLayers, SphLayers - amount of layers in spherical/cylinderical part
%Slices - amount of bezier patches on the circumference in a given
%Layer
%BezierOrder - polynomial order of bezier patches (symmeterical for both
%parameters)

%Output:
%Ncyl - amount of points in a vertical line (defacto z direction) in
%cylinderical part
%Nsph - amount of points in a vertical line (defacto z direction) in
%spherical part
%Ncircm - amount of points in a horizontal (circumference of stump) line

%Ncircum=(BezierOrder+1)+(Slices-1)*BezierOrder=CakeLayers*BezierOrder
%Nvert=(BezierOrder+1)+(Layers-2)*BezierOrder+(BezierOrder-1)=CakeSlices*BezierOrder+1

Ncyl=CylLayers*BezierOrder+1;
Nsph=SphLayers*BezierOrder+1;
Ncircum=Slices*BezierOrder;
end
function CP=StumpMesh2CP(MeshNodes,BezO)
%Input:
%MeshNodes  - point cloud (class double or pointCloud) of size Nvert x Ncircum x 3
%BezO - order of bezier patches (same for both parameters).
%BezO=3 <-> 4 control points in each parameter


%output:
%CP is a struct with:

%CP.v - vertices [x,y,z] ((Nvert*Ncircum)x3) format

%CP.p -  patches. size([BezO+1,BezO+1,Patch Amount]).
%Values are row index of points in CP.v. depth index is patch indexing.

%CP.c - patch connectivity. size([Ptch Amount,4]). The row index of CP.c indicates the patch involved.
%Values are patch indexes. a patch can have 4 neighbors and they are ordered [Left, Right, Up,
%Down] in the row. if a patch has no connectivity in a certin direction,
%the related value will be 0.

%Note1:
%Patches are stiched in the vertical and circumference direction, and there is an overlap of one
%"line" inbetween patches to accomidate g0 geometeric continuity.
% 1st patch,  2nd-(end-1), patch %(end) patch - connects with 1st patch
%   |           |                 |
%   V           V                 V
%
%hence:
%M=Nvert=(BezierOrder+1)+(CakeSlices-2)*BezierOrder+(BezierOrder-1)=CakeSlices*BezierOrder+1
%N=Ncircum=(BezierOrder+1)+(CakeLayers-1)*BezierOrder=CakeLayers*BezierOrder
%where:
%M=Nvert - amount of nodes in a vertical line
%N=Ncircm - amount of nodes in a horizontal (circumference of stump) line

if isa(MeshNodes,'pointCloud'); MeshNodes=MeshNodes.Location; end

sz=size(MeshNodes);
M=sz(1); %amount of nodes in a vertical line
N=sz(2); %amount of nodes in a horizontal (circumference of stump) line
NodeAmnt=M*N;
CakeSlices=(M-1)/BezO; %amount of patches in a single layer
CakeLayers=N/BezO; %amount of layers of patches
PtchAmnt=CakeSlices*CakeLayers;

Iv=reshape(1:NodeAmnt,[M,N]); %indexing matrix to vertices in ControlPoints

%Build V
V=reshape(MeshNodes,[],3); %vertices [x,y,z] (mx3) format

%Build P
P=zeros(BezO+1,BezO+1,PtchAmnt);
k=1;
for n=1:BezO:(N-2*BezO+1) %go "block" by "block" and collect it into P
    for m=1:BezO:(M-BezO)
        P(:,:,k)=Iv(m:m+BezO,n:n+BezO);
        k=k+1;
    end
end
for m=1:BezO:(M-BezO) %add the blocks that stich the circumference togther
    P(:,:,k)=[Iv(m:m+BezO,(N-BezO+1):N),Iv(m:m+BezO,1)];
    k=k+1;
end

%Build C
Ip=reshape(1:PtchAmnt,[CakeLayers,CakeSlices]); %indexing point in PtchCP.
C=zeros(PtchAmnt,4);
for k=1:PtchAmnt
    [i,j]=ind2sub(size(Ip),k); %returns indexes [i,j] for patch number k
    if i==1, U=0; else, U=Ip(i-1,j); end %returns 0 if no patch above
    if i==CakeLayers, D=0; else, D=Ip(i+1,j); end %returns 0 if no patch below
    if j==1, L=Ip(i,CakeSlices); else, L=Ip(i,j-1); end %remember circularity
    if j==CakeSlices, R=Ip(i,1); else, R=Ip(i,j+1); end
    C(k,:)=[L,R,U,D];
    %------------
end

%insert to struct
CP.v=V;
CP.p=P;
CP.c=C;
end
function CP=RadialOptimization4CP(Method,CP,StaticPntCld,Rbounds,Xcenter,varargin)
%Blow up CP.v to so the evaluated bezier surfaces register to StaticPntCld using fmincon (optimization)
%Algorithm works in two steps:
%--rough optimization to fit all CP.v with the same radius of expansion
%--a fine optimization that fits different radii to different vertices

%two methods for radial optimization in this code:
%--cylinderical: radial optimization from the z axis
%--spherical: radial optimization from a given point

%Algorithm Assumption:
%--CP.v are bounded within StaticPntCld.
%--to register beier surfaces onto StaticPntCLd one must expand CP.v radially
%from a given point/axis.

%Example for expansion from a given axis:
%|  <-|   |->  |
%|  <-|   |->  |
%|  <-|   |->  |

%Input:
%Method - 'Cylinderical'/'Spherical'
%StaticPntCld - Point cloud in format of [x,y,z] (mx3) OR pointCloud
%object with similar .Location value
%CP.v - see notes below
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

%Output:
%CP - updated, optimized CP

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
%Default variables
if strcmpi(Method,'Spherical') && (~exist('Xcenter','var') || isempty(Xcenter))
    error('You must include Xcenter if spherical method was chosen. See varargin inputs');
end

%if point clouds give in uncomfi formats, change them
%StaticPntCld needs to be pointCloud object for cost function (use of
%findNearestNeighbors)
if ~isa(StaticPntCld,'pointCloud'), StaticPntCld=pointCloud(StaticPntCld); end


m=size(CP.v,1); %find amount of points to move

%find t - matrix mx3 of radial direction per point
switch lower(Method)
    case 'spherical'
        t=(CP.v-Xcenter)./vecnorm(CP.v-Xcenter,2,2);
    case 'cylinderical'
        xy=CP.v(:,[1,2]);
        t=[xy-Xcenter(1:2),zeros(m,1)]./vecnorm(xy-Xcenter(1:2),2,2);
end

N=size(CP.p,1);
fun=@(r) CostFcn(CP,StaticPntCld,t,r,2*N);
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
CP.v=CP.v+t.*rFine;

    function Cost=CostFcn(CP,StaticPntCld,t,r,N)
        %note: Nested functions variables who share a common name with
        %their parent function will become global. for such reason, 'N' was
        %used instead of 'm'
        Cost=0;
        CP.v=CP.v+t.*r;
        [x,y,z]=CombinePatches(CP,N);
        for i=1:numel(x)
            pnti=[x(i),y(i),z(i)];
            [~,d]=findNearestNeighbors(StaticPntCld,pnti,1);
            Cost=Cost+d^2;
        end
    end
end
%% extra functions (note: can't call one method from another)
function varargout=stlread(file)
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
function [P,Q]=BezBellCurve(R,d,L,Nq)
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
function R=EvalBezCrv_DeCasteljau(Q,q)
%Evaluates Bezier Curve by given nodes and parameter
%value

%Q - nodes in format [x] or [x,y,z] dimension matrix. top row is q=0.
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
function R=EvalBezPtch_DeCasteljau(Pcp,u,v)
%Referance https://pages.mtu.edu/~shene/COURSES/cs3621/NOTES/surface/bezier-de-casteljau.html
%u,v - numeric running values of patch (running 0 to 1)

%Pcp - patch control points will be in the format of size(Pcp)=[Vorder+1,Uorder+1,3]
%for example: control point (1,1) value is placed as such: Pcp(1,1,:) = [X,Y,Z]
%EXAMPLE: Pcp(:,:)= P11,P12,P13 -----> u direction
%                   P21,P22,P23
%                   |
%                   |
%                   |
%     v direction   V

%runs algorithem on control points in the v direction, with parameter v to obtain a set of points.
%then run on those points with paramter u to obtain R(u,v)
Xv=EvalBezCrv_DeCasteljau(Pcp(:,:,1),v);
Yv=EvalBezCrv_DeCasteljau(Pcp(:,:,2),v);
Zv=EvalBezCrv_DeCasteljau(Pcp(:,:,3),v);
R(1)=EvalBezCrv_DeCasteljau(Xv',u);%x
R(2)=EvalBezCrv_DeCasteljau(Yv',u);%y
R(3)=EvalBezCrv_DeCasteljau(Zv',u);%z
end
function [SrfP,U,V]=BezPtchCP2SrfP(Pcp,N)
%retruns evaulted points on bezier patch surface from bezier patch conrol
%points. Amount of points if symmetric for both parametre direction (N)

%N - number of points for patching drawing is N^2
%Pcp- patch control points will be in the format of size(Pcp)=[Vorder+1,Uorder+1,3]
%for example: control point (1,1) value is placed as such: Pcp(1,1,:) = [X,Y,Z]
%EXAMPLE: Pcp(:,:)= P11,P12,P13 -----> u direction
%                   P21,P22,P23
%                   |
%                   |
%                   |
%     v direction   V


%returns Psrfp (Patch surf points) - numeric matrix NxNx3 containing (:,:[x,y,z]) of patch
% and U,V for meshing if anyone wants.
SrfP=zeros(N,N,3); %initalize
u=linspace(0,1,N);
v=linspace(0,1,N);
[U,V]=meshgrid(u,v);
for i=1:N
    for j=1:N
        R=EvalBezPtch_DeCasteljau(Pcp,U(i,j),V(i,j)); %Obtain Point
        SrfP(i,j,1)=R(1); SrfP(i,j,2)=R(2); SrfP(i,j,3)=R(3); %Place in Matrix
    end
end
end
function [z,r,residual]=fitcircle(x,varargin)
% Obtained from "fitcircle.m" from mathworks file exchange
% https://www.mathworks.com/matlabcentral/fileexchange/15060-fitcircle-m

%FITCIRCLE    least squares circle fit
%
%   [Z, R] = FITCIRCLE(X) fits a circle to the N points in X minimising
%   geometric error (sum of squared distances from the points to the fitted
%   circle) using nonlinear least squares (Gauss Newton)
%       Input
%           X : 2xN array of N 2D points, with N >= 3
%       Output
%           Z : center of the fitted circle
%           R : radius of the fitted circle
%
%   [Z, R] = FITCIRCLE(X, 'linear') fits a circle using linear least
%   squares minimising the algebraic error (residual from fitting system
%   of the form ax'x + b'x + c = 0)
%
%   [Z, R] = FITCIRCLE(X, Property, Value, ...) allows parameters to be
%   passed to the internal Gauss Newton method. Property names can be
%   supplied as any unambiguous contraction of the property name and are
%   case insensitive, e.g. FITCIRCLE(X, 't', 1e-4) is equivalent to
%   FITCIRCLE(X, 'tol', 1e-4). Valid properties are:
%
%       Property:                 Value:
%       --------------------------------
%       maxits                    positive integer, default 100
%           Sets the maximum number of iterations of the Gauss Newton
%           method
%
%       tol                       positive constant, default 1e-5
%           Gauss Newton converges when the relative change in the solution
%           is less than tol
%
%   [X, R, RES] = fitcircle(...) returns the 2 norm of the residual from
%   the least squares fit
%
%   Example:
%       x = [1 2 5 7 9 3; 7 6 8 7 5 7];
%       % Get linear least squares fit
%       [zl, rl] = fitcircle(x, 'linear')
%       % Get true best fit
%       [z, r] = fitcircle(x)
%
%   Reference: "Least-squares fitting of circles and ellipses", W. Gander,
%   G. Golub, R. Strebel - BIT Numerical Mathematics, 1994, Springer

% This implementation copyright Richard Brown, 2007, but is freely
% available to copy, use, or modify as long as this line is maintained

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
function [K,H,Pmax,Pmin]=SurfCurvature(X,Y,Z)
%Created by Daniel Claxton
%Taken from 
% https://www.mathworks.com/matlabcentral/fileexchange/11168-surface-curvature

% SURFATURE -  COMPUTE GAUSSIAN AND MEAN CURVATURES OF A SURFACE
%   [K,H] = SURFATURE(X,Y,Z), WHERE X,Y,Z ARE 2D ARRAYS OF POINTS ON THE
%   SURFACE.  K AND H ARE THE GAUSSIAN AND MEAN CURVATURES, RESPECTIVELY.
%   SURFATURE RETURNS 2 ADDITIONAL ARGUEMENTS,
%   [K,H,Pmax,Pmin] = SURFATURE(...), WHERE Pmax AND Pmin ARE THE MINIMUM
%   AND MAXIMUM CURVATURES AT EACH POINT, RESPECTIVELY.


% Example: 
% [X,Y,Z] = peaks; 
% [K,H,P1,P2] = SurfCurvature(X,Y,Z); 
% surf(X,Y,Z,H,'facecolor','interp'); 
% set(gca,'clim',[-1,1])

% First Derivatives
[Xu,Xv] = gradient(X);
[Yu,Yv] = gradient(Y);
[Zu,Zv] = gradient(Z);

% Second Derivatives
[Xuu,Xuv] = gradient(Xu);
[Yuu,Yuv] = gradient(Yu);
[Zuu,Zuv] = gradient(Zu);

[Xuv,Xvv] = gradient(Xv);
[Yuv,Yvv] = gradient(Yv);
[Zuv,Zvv] = gradient(Zv);

% Reshape 2D Arrays into Vectors
Xu = Xu(:);   Yu = Yu(:);   Zu = Zu(:); 
Xv = Xv(:);   Yv = Yv(:);   Zv = Zv(:); 
Xuu = Xuu(:); Yuu = Yuu(:); Zuu = Zuu(:); 
Xuv = Xuv(:); Yuv = Yuv(:); Zuv = Zuv(:); 
Xvv = Xvv(:); Yvv = Yvv(:); Zvv = Zvv(:); 

Xu          =   [Xu Yu Zu];
Xv          =   [Xv Yv Zv];
Xuu         =   [Xuu Yuu Zuu];
Xuv         =   [Xuv Yuv Zuv];
Xvv         =   [Xvv Yvv Zvv];

% First fundamental Coeffecients of the surface (E,F,G)
E           =   dot(Xu,Xu,2);
F           =   dot(Xu,Xv,2);
G           =   dot(Xv,Xv,2);

m           =   cross(Xu,Xv,2);
p           =   sqrt(dot(m,m,2));
n           =   m./[p p p]; 

% Second fundamental Coeffecients of the surface (L,M,N)
L           =   dot(Xuu,n,2);
M           =   dot(Xuv,n,2);
N           =   dot(Xvv,n,2);

[s,t] = size(Z);

% Gaussian Curvature
K = (L.*N - M.^2)./(E.*G - F.^2);
K = reshape(K,s,t);

% Mean Curvature
H = (E.*N + G.*L - 2.*F.*M)./(2*(E.*G - F.^2));
H = reshape(H,s,t);

% Principal Curvatures
Pmax = H + sqrt(H.^2 - K);
Pmin = H - sqrt(H.^2 - K);
end
function [x,y,z]=CombinePatches(CP,N)
%Input:
%CP is a struct with:

%CP.v - vertices [x,y,z] ((Nvert*Ncircum)x3) format

%CP.p -  patches. size([BezO+1,BezO+1,Patch Amount]).
%Values are row index of points in CP.v. depth index is patch indexing.

%CP.c - patch connectivity. size([Ptch Amount,4]). The row index of CP.c indicates the patch involved.
%Values are patch indexes. a patch can have 4 neighbors and they are ordered [Left, Right, Up,
%Down] in the row. if a patch has no connectivity in a certin direction,
%the related value will be 0.

%N - number of points evaluated for each patch is N^2

%Output:
%x,y,z matrices size PtchAmnt*N x N of all patches combined


%surface produced by [x,y,z] is distorted. no real advantage

%Input:
%CP struct containing v,p,c matrices
%N - number of points evaluated per patch are N^2

%output:
%x,y,z matrices containing data of all patches combined.
PtchAmnt=size(CP.p,3);
sz=size(CP.p);
[x,y,z]=deal(zeros(PtchAmnt*N,N));
for k=1:PtchAmnt
    Pcp=reshape(CP.v(CP.p(:,:,k),:),[sz(1:2),3]); %Obtain patch control points. format explained above
    Psrfp=BezPtchCP2SrfP(Pcp,N);
    x((k-1)*N+1:k*N,:,:)=Psrfp(:,:,1);
    y((k-1)*N+1:k*N,:,:)=Psrfp(:,:,2);
    z((k-1)*N+1:k*N,:,:)=Psrfp(:,:,3);
end
end
%% functions not in use - ideas not implemented
function [xyz,Psi,Theta]=CalcSphere(R,Npsi,Ntheta,varargin)
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
%INPUT:
%R - raidus
%Ntheta - amount of theta points (descrization). theta ranges [0,2pi]
%Nz -  amount of z points (descrization). z ranges [0,L]

%Varagin:
%Trans - translation [x,y,z]

%OUTPUT:
%xyz - size Nz x Ntheta x 3 where the depth dimension is constructed as [x,y,z]
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
function r=FindSphereRadius(xyz)
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
zy=xyz(:,[2,3]);
zyhull=zy(convhull(zy),:);
[~,rzy]=fitcircle(zyhull,'linear');

zx=xyz(:,[1,3]);
zxhull=zy(convhull(zx),:);
[~,rzx]=fitcircle(zxhull,'linear');

r=mean([rzx,rzy]);
end
function BlownPntCld=RadialNonRigidRegistration(Method,MovPntCld,StaticPntCld,Rbounds,Xcenter,varargin)
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