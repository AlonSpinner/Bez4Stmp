classdef BezCP
    properties (SetAccess = private)
        BezierOrder
        %integer. symmetric for both u and v parameters
        Patches
        %integers size BezierOrder+1 x BezierOrder+1 x PatchAmount.
        %Values are row index of points in CP.v. depth index is patch indexing.
        Connectivity
        %Patch connectivity. size([Ptch Amount,4]). The row index of CP.c indicates the patch involved.
        %Values are patch indexes. a patch can have 4 neighbors and they are ordered [Left, Right, Up,
        %Down] in the row. if a patch has no connectivity in a certin direction,
        %the related value will be 0.
    end
    properties
        Vertices
        %double size mx3 of control points
    end
    methods %public methods, accessible from outside the class
        function obj=BezCP(MeshNodes,BezO,varargin) %constructor
            %Input:
            %MeshNodes  - point cloud (class double or pointCloud) of size N x M x 3
            %BezO - order of bezier patches (same for both parameters).
            %BezO=3 <-> 4 control points in each parameter
            
            %Varargin Input:
            %Method - 'Circular'/'CircularWithCap'/'Block'. Block
            %Circular - stich patches to close a cylinderical surface.
            %stiching occures across columns of MeshNodes:
            %first patch and last patch both consist of
            %MeshNodes(:,1)
            %M=Nvert=(BezierOrder+1)+(Layers-2)*BezierOrder+(BezierOrder-1)=Layers*BezierOrder+1
            %N=Ncircum=(BezierOrder+1)+(Slices-1)*BezierOrder=Slices*BezierOrder
            %Block
            %M=Nvert=Layers*BezierOrder+1
            %N=Nhorz=Slices*BezierOrder+1
            
            %Note:
            %Patches are stiched in the vertical and horizontal direction, and there is an overlap of one
            %"line" inbetween patches to accomidate g0 geometeric continuity.
            % 1st patch,  2nd-(end-1), patch %(end) patch - connects with 1st patch
            %   |           |                 |
            %   V           V                 V
            %
            %M=Nvert - amount of nodes in a vertical line
            %N=Nhorz - amount of nodes in a horizontal line
            
            if isa(MeshNodes,'pointCloud'); MeshNodes=MeshNodes.Location; end
            
            Method='block'; %default
            for ind=1:2:length(varargin)
                comm=lower(varargin{ind});
                switch comm
                    case 'method'
                        Method=varargin{ind+1};
                end
            end
            
            sz=size(MeshNodes);
            M=sz(1); %amount of nodes in a vertical line
            N=sz(2); %amount of nodes in a horizontal line
            NodeAmnt=M*N;
            
            %Build V - vertices
            Iv=reshape(1:NodeAmnt,[M,N]); %indexing matrix to vertices in ControlPoints
            V=reshape(MeshNodes,[],3); %vertices [x,y,z] (mx3) format
            
            switch lower(Method)
                case 'circularwithcap'
                    Layers=(M-1)/BezO; %amount of patches in a single layer
                    Slices=N/BezO; %amount of layers of patches
                    %check validity
                    if mod(Layers,1) || mod(Slices,1)
                        error(['for the given method, (M-1)/BezO and N/BezO must both be integers',...
                            'where MeshNodes is of size N x M x 3']);
                    end
                    if 4~=Slices
                        error('for the given method, Slices must equal 4')
                    end
                    PtchAmnt=Slices*Layers+1; %+1 for cap
                    
                    %add cap vertices to V
                    CapCircum=Iv(end,:); %indcies of cap circumference nodes in V
                    Vcap=zeros((BezO-1)^2,3); %initalize vertices of cap
                    for j=2:BezO %linear interpolation
                        for i=2:BezO
                            p1=V(CapCircum(BezO+j),:);
                            p2=V(CapCircum(end-j+2),:);
                            Vcap((i-1)*(j-1),:)=(p1*i+p2*(BezO+1-i))/(BezO+1);
                        end
                    end
                    V=[V;Vcap]; %add to V
                    
                    %Build P - patches
                    P=zeros(BezO+1,BezO+1,PtchAmnt);
                    k=1;
                    for n=1:BezO:(N-2*BezO+1) %go "block" by "block" and collect it into P
                        for m=1:BezO:(M-BezO)
                            P(:,:,k)=Iv(m:m+BezO,n:n+BezO);
                            k=k+1;
                        end
                    end
                    
                    %add the blocks that stich the circumference togther
                    for m=1:BezO:(M-BezO)
                        P(:,:,k)=[Iv(m:m+BezO,(N-BezO+1):N),Iv(m:m+BezO,1)];
                        k=k+1;
                    end
                    
                    %add cap
                    %sides
                    P(1:BezO+1,1,k)=CapCircum(1:BezO+1); %left
                    P(end,2:end-1,k)=CapCircum(BezO+2:2*BezO); %buttom
                    P(1:BezO+1,end,k)=fliplr(CapCircum(2*BezO+1:3*BezO+1)); %right
                    P(1,2:end-1,k)=fliplr(CapCircum(3*BezO+2:end)); %top
                    %middle
                    for j=2:BezO
                        for i=2:BezO
                            P(i,j,k)=(i-1)*(j-1)+NodeAmnt;
                        end
                    end
                    
                    %Build C - connectivity
                    %cap patch has no row in C, but is indexed there
                    Ip=reshape(1:PtchAmnt-1,[Layers,Slices]); %indexing point in PtchCP. not including cap
                    Icap=PtchAmnt; %index number of cap patch
                    C=zeros(PtchAmnt-1,4);
                    for k=1:(PtchAmnt-1) %run on all patches not including cap
                        [i,j]=ind2sub(size(Ip),k); %returns indexes [i,j] for patch number k
                        if i==1, U=Icap; else, U=Ip(i-1,j); end %returns cap patch index if @ top edge
                        if i==Layers, D=0; else, D=Ip(i+1,j); end %returns 0 if no patch below
                        if j==1, L=Ip(i,Slices); else, L=Ip(i,j-1); end %remember circularity
                        if j==Slices, R=Ip(i,1); else, R=Ip(i,j+1); end %remember circularity
                        C(k,:)=[L,R,U,D];
                    end
                    
                case 'circular'
                    Layers=(M-1)/BezO; %amount of patches in a single layer
                    Slices=N/BezO; %amount of layers of patches
                    %check validity
                    if mod(Layers,1) || mod(Slices,1)
                        error(['for the given method, (M-1)/BezO and N/BezO must both be integers',...
                            'where MeshNodes is of size N x M x 3']);
                    end
                    PtchAmnt=Slices*Layers;
                    
                    %Build P - patches
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
                    
                    %Build C - connectivity
                    Ip=reshape(1:PtchAmnt,[Layers,Slices]); %indexing point in PtchCP.
                    C=zeros(PtchAmnt,4);
                    for k=1:PtchAmnt
                        [i,j]=ind2sub(size(Ip),k); %returns indexes [i,j] for patch number k
                        if i==1, U=0; else, U=Ip(i-1,j); end %returns 0 if no patch above
                        if i==Layers, D=0; else, D=Ip(i+1,j); end %returns 0 if no patch below
                        if j==1, L=Ip(i,Slices); else, L=Ip(i,j-1); end %remember circularity
                        if j==Slices, R=Ip(i,1); else, R=Ip(i,j+1); end %remember circularity
                        C(k,:)=[L,R,U,D];
                    end
                    
                case 'block'
                    Layers=(M-1)/BezO;
                    Slices=(N-1)/BezO;
                    %check validity
                    if mod(Layers,1) || mod(Slices,1)
                        error(['for the given method, (M-1)/BezO and (N-1)/BezO must both be integers',...
                            'where MeshNodes is of size N x M x 3']);
                    end
                    PtchAmnt=Slices*Layers;
                    
                    %Build P - patches
                    P=zeros(BezO+1,BezO+1,PtchAmnt);
                    k=1;
                    for n=1:BezO:(N-2*BezO+1) %go "block" by "block" and collect it into P
                        for m=1:BezO:(M-BezO)
                            P(:,:,k)=Iv(m:m+BezO,n:n+BezO);
                            k=k+1;
                        end
                    end
                    
                    %Build C - connectivity
                    Ip=reshape(1:PtchAmnt,[Layers,Slices]); %indexing point in PtchCP.
                    C=zeros(PtchAmnt,4);
                    for k=1:PtchAmnt
                        [i,j]=ind2sub(size(Ip),k); %returns indexes [i,j] for patch number k
                        if i==1, U=0; else, U=Ip(i-1,j); end %returns 0 if no patch above
                        if i==Layers, D=0; else, D=Ip(i+1,j); end %returns 0 if no patch below
                        if j==1, L=Ip(i,Slices); else, L=0; end %return 0 if no patch to the left
                        if j==Slices, R=Ip(i,1); else, R=0; end  %return 0 if no patch to the right
                        C(k,:)=[L,R,U,D];
                    end
            end
            
            %Introduce to object
            obj.BezierOrder=BezO;
            obj.Vertices=V;
            obj.Patches=P;
            obj.Connectivity=C;
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
            PtchAmnt=size(obj.Patches,3);
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
            %Draw patches
            for k=1:PtchAmnt
                %Evaluate points of patch
                Pcp=reshape(obj.Vertices(obj.Patches(:,:,k),:),[obj.BezierOrder+1,obj.BezierOrder+1,3]); %patch control points.
                Psrfp=BezPtchCP2SrfP(Pcp,N); %compute surface points
                x=Psrfp(:,:,1); y=Psrfp(:,:,2); z=Psrfp(:,:,3);
                
                %Draw evaluated points with surf command
                %Note: switch coded in forloop to allow cool "build up" graphics when plotting
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
                        H=filloutliers(H,'nearest');
                        Handles(k)=surf(Ax,x,y,z,H,'facecolor','interp','edgecolor','none');
                end
                pause(PauseTime);
            end
        end
        function xyz=CombinePatches(obj,N)
            %N - number of points evaluated for each patch is N^2
            
            %Output:
            %xyz size PtchAmnt*N x N x 3 of all patches combined where the
            %third dimension is ordered [x,y,z]
            
            PtchAmnt=size(obj.Patches,3);
            sz=size(obj.Patches);
            [x,y,z]=deal(zeros(PtchAmnt*N,N));
            for k=1:PtchAmnt
                Pcp=obj.Vertices(obj.Patches(:,:,k),:); %Obtain patch control points.
                Pcp=reshape(Pcp,[sz(1:2),3]); %format to match requirements of BezPtchCP2SrfP
                Psrfp=BezPtchCP2SrfP(Pcp,N);
                x((k-1)*N+1:k*N,:,:)=Psrfp(:,:,1);
                y((k-1)*N+1:k*N,:,:)=Psrfp(:,:,2);
                z((k-1)*N+1:k*N,:,:)=Psrfp(:,:,3);
            end
            xyz=cat(3,x,y,z);
        end
        
        function igsWrite(obj,FileName)
            %Input:
            %FileName - char. can contain path. may not include extension.
            
            %Output:
            %creates iges file format by the following steps:
                %1) creates itd text format and writes control point data
                    %into it
                %2) call on irit2igs.exe to create igs file from itd
                %3) delete itd file
                
            %NOTE:
            %can only run if file 'irit2igs.exe' is in folder/matlab path
                
            
            %example of itd text format:
            %[SURFACE BEZIER 3 3 E3
            % 	[-42.6193	-3.02585	1.23967	] %control points for u=const1
            % 	[-42.2107	-3.74057	5.48943	]
            % 	[-41.9059	-3.41341	9.40101	]
            %
            % 	[-42.1267	-10.0367	1.25626	] %control points for u=const2
            % 	[-41.6183	-9.73169	5.64561	]
            % 	[-41.1738	-10.0765	9.31034	]
            %
            % 	[-40.3078	-16.7856	1.19208	] %control points for u=const3
            % 	[-39.9241	-16.3549	5.55235	]
            % 	[-39.1993	-16.8461	9.42853	]
            %]
            
            %create FileNames for igs and itd
            %enforce '.itd' extention to FileName if other was provided
            [~,~,ext]=fileparts(FileName);
            FileName_itd=[FileName(1:end-length(ext)) '.itd'];
            FileName_igs=[FileName(1:end-length(ext)) '.igs'];
            
            %open itd file to write to
            fileID = fopen(FileName_itd,'w');
            if fileID<0 %=-1
                error('FileName is invalid');
            end
            
            %initialize. m=n but seprated for clarity
            PatchAmnt=size(obj.Patches,3);
            m=size(obj.Patches,1);
            n=size(obj.Patches,2);
            
            for k1=1:PatchAmnt %amount of patches
                P=obj.Vertices(obj.Patches(:,:,k1),:); %collect patch points in m*nx3 matrix
                fprintf(fileID,'[SURFACE BEZIER %d %d E3\n',m,n); %open line for patch k1
                for k2=1:m:(m*n-m+1) %every m rows represent a "block" (constant u)
                    for k3=0:(n-1) %print block
                        fprintf(fileID,['\t[',sprintf('%g\t',P(k2+k3,:)),']\n']);
                    end
                    if k2~=(m*n-m+1), fprintf(fileID,'\n'); end %go down a line between blocks
                end
                fprintf(fileID,']\n\n'); %close line for patch k1
            end
            fclose(fileID); %close connection between file and matlab
            
            %use irit to convert itd file to igs
            eval(sprintf('!irit2igs -o %s %s',FileName_igs,FileName_itd));
            
            %delete itd file
            delete(FileName_itd);
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
    end
end
%% extra functions (note: can't call one method from another)
%Bezier evaluating functions
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
%curvature evaluation
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