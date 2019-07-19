%{
Written by Alon Spinner @ 6-7/19
handles.GUIstruct is a struct with fieldnames
-Scan: mx3 of class double describing a point cloud
-Stmp: Bez4Stmp class object
%}
function varargout = Bez4StmpGUI(varargin)
% Bez4StmpGUI MATLAB code for Bez4StmpGUI.fig
%      Bez4StmpGUI, by itself, creates a new Bez4StmpGUI or raises the existing
%      singleton*.
%
%      H = Bez4StmpGUI returns the handle to a new Bez4StmpGUI or the handle to
%      the existing singleton*.
%
%      Bez4StmpGUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in Bez4StmpGUI.M with the given input arguments.
%
%      Bez4StmpGUI('Property','Value',...) creates a new Bez4StmpGUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before Bez4StmpGUI_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to Bez4StmpGUI_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help Bez4StmpGUI

% Last Modified by GUIDE v2.5 18-Jul-2019 18:51:15

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
    'gui_Singleton',  gui_Singleton, ...
    'gui_OpeningFcn', @Bez4StmpGUI_OpeningFcn, ...
    'gui_OutputFcn',  @Bez4StmpGUI_OutputFcn, ...
    'gui_LayoutFcn',  [] , ...
    'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
function Bez4StmpGUI_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to Bez4StmpGUI (see VARARGIN)

% Choose default command line output for Bez4StmpGUI
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

%establish Ax to draw on
Ax=handles.Ax;
Ax=BezCP.CreateDrawingAxes(Ax);

%Create BezCP to draw with
PtchAmnt=16;
Layers=sqrt(PtchAmnt);
BezO=3;
N=Layers*BezO+1;
[X,Y,Z]=peaks(N);
X=X*3; Y=Y*3; %scale
MeshNodes=zeros(size(X));
MeshNodes(:,:,1)=X; MeshNodes(:,:,2)=Y; MeshNodes(:,:,3)=Z;
CP=BezCP(MeshNodes,BezO,'Method','Block');
PseudoInverseCP=CP.PesudoInverseVertices;

%Draw
PseudoInverseCP.DrawBezierPatches('Ax',Ax,'facealpha',1,...
    'edgecolor',0.5*[1,1,1],'N',10);
BezCP.DrawPointCloud(PseudoInverseCP.Vertices,'Ax',Ax,'color',[0,1,0],'msize',20); %draw control points
BezCP.DrawPointCloud(MeshNodes,'Ax',Ax,'color',[1,0,0],'msize',20);

%turn rotation off. For some reason "pcshow" function (called in
%"BezCP.DrawPointCloud" turns rotate3d on for axes instilled in guide GUI.
rotate3d(Ax,'off');
handles.RotateToolToggle.State='off';
function varargout = Bez4StmpGUI_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure

varargout{1} = handles.output;
%% Callbacks
%action related
function FileVarNameEdit_Callback(hObject, eventdata, handles)
%cancel all actions but load (turn pushbuttons dark)
DarkGreen=[0.31,0.49,0.37]; DarkBlue=[0.4,0.58,0.62];
handles.ScanPointCloudPush.BackgroundColor=DarkGreen;
handles.BezierSurfaceMeshPlotPush.BackgroundColor=DarkGreen;
handles.HausdorffPlotPush.BackgroundColor=DarkGreen;
handles.CalculatePush.BackgroundColor=DarkBlue;
handles.ExportPush.BackgroundColor=DarkBlue;
handles.SaveWorkspacePush.BackgroundColor=DarkBlue;
function LoadPush_Callback(hObject, eventdata, handles)
%{
Updates handles.GUIstruct.Scan by user input and drawns it unto handles.Ax

if file chosen in FileVarNamePopUp:
accepts .stl and .m files
.m files should contain Bez4Stmp/pointCloud/double class having the same name
as the file they were in. if they dont, chooses by hierarchy
(Bez4Stmp>pointCloud>double)

if varaible chosen in FileVarNamePopUp:
looks for Bez4Stmp/pointCloud/double class with the given name in workspace
and loads it.
%}

%clear axes incase of error
Ax=handles.Ax;
cla(Ax,'reset');
BezCP.CreateDrawingAxes(Ax);

str=handles.FileVarNamePopUp.String{handles.FileVarNamePopUp.Value};
switch str
    case 'File with extension'
        FileFullName=handles.FileVarNameEdit.String;
        [~,FileName,ext]=fileparts(FileFullName);
        switch ext
            case '.stl'
                [~,Scan]=stlread(FileFullName); %obtain vertices of STL
            case '.mat'
                S=load(FileFullName); %S is a struct with fieldnames being the names of the variables in mat
                NameMatches=strcmp(fieldnames(S),FileName); %logical array
                cellS=struct2cell(S); %cell array containig struct fields (loses field names)
                Bez4StmpClassMatches=cellfun(@(s) isa(s,'Bez4Stmp'),cellS); %logical array
                doubleClassMatches=cellfun(@(s) isa(s,'double'),cellS); %logical array
                pointCloudClassMatches=cellfun(@(s) isa(s,'pointCloud'),cellS); %logical array
                if any(NameMatches) %check for a name match (maximum 1 possible)
                    var=cellS{NameMatches};
                    switch class(var)
                        case 'Bez4Stmp'
                            Stmp=var;
                            Scan=Stmp.Scan.Location;
                        case 'pointCloud'
                            Scan=var.Location;
                        case 'double'
                            Scan=var;
                    end
                    %if no variable with file name exists, take first variable by
                    %preference: Bez4Stmp>pointCloud>double
                elseif any(Bez4StmpClassMatches) %Bez4Stmp
                    Stmp=cellS{Bez4StmpClassMatches(1)};
                    Scan=Stmp.Scan;
                elseif any(pointCloudClassMatches) %pointCloud
                    Scan=cellS{pointCloudClassMatches(1)}.Location;
                elseif any(doubleClassMatches) %double
                    Scan=cellS{doubleClassMatches(1)};
                else %.m file contains no good data
                    error(sprintf(['Wrong Input\n',...
                        'Input file must be a .m or .stl one\n',...
                        '.m files will contain a Bez4Stmp, double or pointCloud class object.\n',...
                        'It is recommended the the variable name have the same name as the .m file']));
                end %end .m related if-elseif
        end %end .ext related switch-case
    case 'Workspace variable'
        VarName=handles.FileVarNameEdit.String;
        Var=evalin('base',VarName);
        switch class(Var)
            case 'pointCloud'
                Scan=pointCloud.Location;
            case 'double'
                Scan=Var;
            case 'Bez4Stmp'
                Stmp=Var;
                Scan=Stmp.Scan.Location;
        end
end %switch between workspace variable / file

%do some further input checking for pointCloud/double inputs
sz=size(Scan);
if numel(sz)~=2 || sz(2)~=3
    error('Input pointCloud.Location or double matrix must have size mx3');
end

%enable/disable specific actions (if color is DarkGreen/DarkBlue buttons are disabled and otherwise)
%update GUIstruct Scan and Stmp (if Bez4Stmp was loaded)
LightGreen=[0.43,0.8,0.55]; DarkGreen=[0.31,0.49,0.37];
LightBlue=[0.53,0.8,0.86]; DarkBlue=[0.4,0.58,0.62];
if exist('Stmp','var') %if variable "Stmp" was created (Bez4Stmp input)
    handles.GUIstruct.Stmp=Stmp; %insert Stmp into GUIstruct to be called by other callbacks
    handles.GUIstruct.Scan=Scan; %insert Scan into GUIstruct to be called by other callbacks
    
    handles.ScanPointCloudPush.BackgroundColor=LightGreen;
    handles.BezierSurfaceMeshPlotPush.BackgroundColor=LightGreen;
    handles.HausdorffPlotPush.BackgroundColor=LightGreen;
    handles.CurvaturePlotPush.BackgroundColor=LightGreen;
    handles.CalculatePush.BackgroundColor=LightBlue;
    handles.ExportPush.BackgroundColor=LightBlue;
    handles.SaveWorkspacePush.BackgroundColor=LightBlue;
else %if variable "Stmp" wasnt created (double or pointCloud input)
    handles.GUIstruct.Scan=Scan; %insert Scan into GUIstruct to be called by other callbacks
    
    handles.ScanPointCloudPush.BackgroundColor=LightGreen;
    handles.CalculatePush.BackgroundColor=LightBlue;
    
    handles.BezierSurfaceMeshPlotPush.BackgroundColor=DarkGreen;
    handles.HausdorffPlotPush.BackgroundColor=DarkGreen;
    handles.CurvaturePlotPush.BackgroundColor=DarkGreen;
    handles.ExportPush.BackgroundColor=DarkBlue;
    handles.SaveWorkspacePush.BackgroundColor=DarkBlue;
end

guidata(hObject, handles); %save handles in figure
ScanPointCloudPush_Callback(handles.ScanPointCloudPush,[],handles)%Draw InputPointCloud
function CalculatePush_Callback(hObject, eventdata, handles)
DarkBlue=[0.4,0.58,0.62];
if norm(hObject.BackgroundColor-DarkBlue)<0.01
    errordlg('Please load data prior');
    return
end

handles=guidata(hObject); %Obtain updated handles

%Obtain data from Edit
CylLayers=str2double(handles.CylinderLayersEdit.String);
SphLayers=str2double(handles.SphereLayersEdit.String);
Slices=str2double(handles.SlicesEdit.String);
BezierOrder=str2double(handles.BezierOrderEdit.String);
Cap=str2num(handles.CapPopUp.String{handles.CapPopUp.Value});
XcenterCalcMethod=handles.XcenterPopUp.String{handles.XcenterPopUp.Value};

%compute stmp
Scan=handles.GUIstruct.Scan;
Stmp=Bez4Stmp(Scan,'Cap',Cap,'SphLayers',SphLayers,'CylLayers',CylLayers,...
    'Slices',Slices,'BezierOrder',BezierOrder,'XcenterCalculationMethod',XcenterCalcMethod);

%alert user of succsessful load and enable actions by changing pushbuttons
%color
helpdlg('Calculation process has been completed');
LightGreen=[0.43,0.8,0.55]; LightBlue=[0.53,0.8,0.86];
handles.BezierSurfaceMeshPlotPush.BackgroundColor=LightGreen;
handles.HausdorffPlotPush.BackgroundColor=LightGreen;
handles.CurvaturePlotPush.BackgroundColor=LightGreen;
handles.ExportPush.BackgroundColor=LightBlue;
handles.SaveWorkspacePush.BackgroundColor=LightBlue;

handles.GUIstruct.Stmp=Stmp; %save to GUI struct
guidata(hObject, handles); %save handles in figure
function SaveWorkspacePush_Callback(hObject, eventdata, handles)
DarkBlue=[0.4,0.58,0.62];
if norm(hObject.BackgroundColor-DarkBlue)<0.01
    errordlg('Please calculate data prior');
    return
end

handles=guidata(hObject); %Obtain updated handles
Stmp=handles.GUIstruct.Stmp;
Uinput=inputdlg('Please enter variable name','',1,{'MyStump'});
StmpName=Uinput{1};
assignin('base',StmpName,Stmp);
function ExportPush_Callback(hObject, eventdata, handles)
DarkBlue=[0.4,0.58,0.62];
if norm(hObject.BackgroundColor-DarkBlue)<0.01
    errordlg('Please calculate data prior');
    return
end

handles=guidata(hObject); %Obtain updated handles
Stmp=handles.GUIstruct.Stmp;
filter={'*.igs';'*.mat'};
[file,path]=uiputfile(filter); %file=[name,ext]
if ~ischar(file), return, end %user cancelled operation in ui
[~,name,ext]=fileparts(file);
switch ext
    case '.igs'
        Stmp.StmpBezCP.igsWrite(name);
    case '.mat' %"dynamic fieldname technique"
        S.(name)=Stmp; %create struct with field [name chosen by user]
        save([path,file],'-struct','S') %save .mat file with variable named [name chosen by user]
end
helpdlg(sprintf('%s file saved succesfully',ext));
%plotting callbacks
function HausdorffPlotPush_Callback(hObject, eventdata, handles)
DarkGreen=[0.31,0.49,0.37];
if norm(hObject.BackgroundColor-DarkGreen)<0.01
    errordlg('Please calculate data prior');
    return
end

handles=guidata(hObject); %Obtain updated handles
Stmp=handles.GUIstruct.Stmp;

%restablish Ax to draw on
Ax=handles.Ax;
cla(Ax,'reset');
Ax=BezCP.CreateDrawingAxes(Ax);

%Draw
Zfilter=str2num(handles.HausedorffZThresholdEdit.String);
Stmp.HausdorffAsses('Ax',Ax,'zthreshold',Zfilter);

%turn rotation off. For some reason "pcshow" function (called in
%"BezCP.DrawPointCloud" turns rotate3d on for axes instilled in guide GUI.
rotate3d(Ax, 'off');
function BezierSurfaceMeshPlotPush_Callback(hObject, eventdata, handles)
DarkGreen=[0.31,0.49,0.37];
if norm(hObject.BackgroundColor-DarkGreen)<0.01
    errordlg('Please calculate data prior');
    return
end

handles=guidata(hObject); %Obtain updated handles
Stmp=handles.GUIstruct.Stmp;

%restablish Ax to draw on
Ax=handles.Ax;
cla(Ax,'reset');
Ax=BezCP.CreateDrawingAxes(Ax);

%Draw
Stmp.StmpBezCP.DrawBezierPatches('Ax',Ax);
function ScanPointCloudPush_Callback(hObject, eventdata, handles)
DarkGreen=[0.31,0.49,0.37];
if norm(hObject.BackgroundColor-DarkGreen)<0.01
    errordlg('Please calculate data prior');
    return
end

handles=guidata(hObject); %Obtain updated handles
Scan=handles.GUIstruct.Scan;

%restablish Ax to draw on
Ax=handles.Ax;
cla(Ax,'reset');
Ax=BezCP.CreateDrawingAxes(Ax);

%Draw
BezCP.DrawPointCloud(Scan,'Ax',Ax);

%turn rotation off. For some reason "pcshow" function (called in
%"BezCP.DrawPointCloud" turns rotate3d on for axes instilled in guide GUI.
rotate3d(Ax, 'off');
function CurvaturePlotPush_Callback(hObject, eventdata, handles)
DarkGreen=[0.31,0.49,0.37];
if norm(hObject.BackgroundColor-DarkGreen)<0.01
    errordlg('Please calculate data prior');
    return
end

handles=guidata(hObject); %Obtain updated handles
Stmp=handles.GUIstruct.Stmp;

%restablish Ax to draw on
Ax=handles.Ax;
cla(Ax,'reset');
Ax=BezCP.CreateDrawingAxes(Ax);
colormap(Ax,'jet');k
colorbar(Ax);

%Draw. Currently discrete gives better results for some reason
Type=handles.CurvaturePopUp.String{handles.CurvaturePopUp.Value};
Stmp.StmpBezCP.DrawMeshCurvature('Type',Type,'Ax',Ax,'Technique','Discrete');

%info callbacks
function DocumentationPush_Callback(hObject, eventdata, handles)
%% Functions
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
%% GUI functions with little use (prehaps only input check)
function CylinderLayersEdit_Callback(hObject, eventdata, handles)
[Val,status]=str2num(hObject.String);
if ~status || mod(Val,1)~=0 || ~(Val>0), errordlg('Input Cylinder Layers must be an natrual number'); end
function SphereLayersEdit_Callback(hObject, eventdata, handles)
[Val,status]=str2num(hObject.String);
if ~status || mod(Val,1)~=0 || ~(Val>0), errordlg('Input Sphere Layers must be an natrual number'); end
function SlicesEdit_Callback(hObject, eventdata, handles)
[Val,status]=str2num(hObject.String);
if ~status || mod(Val,1)~=0 || ~(Val>0), errordlg('Input Slices must be an natrual number'); end
function BezierOrderEdit_Callback(hObject, eventdata, handles)
[Val,status]=str2num(hObject.String);
if ~status || mod(Val,1)~=0 || ~(Val>0), errordlg('Input BezierOrder must be an natrual number'); end
function HausedorffZThresholdEdit_Callback(hObject, eventdata, handles)
[Val,status]=str2num(hObject.String);
if ~status, errordlg('Input must be numeric'); end
%% GUI functions with no use
function CapPopUp_Callback(hObject, eventdata, handles)
function CapPopUp_CreateFcn(hObject, eventdata, handles)
% hObject    handle to CapPopUp (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function XcenterPopUp_Callback(hObject, eventdata, handles)
function XcenterPopUp_CreateFcn(hObject, eventdata, handles)
% hObject    handle to XcenterPopUp (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function FileVarNamePopUp_Callback(hObject, eventdata, handles)
function FileVarNamePopUp_CreateFcn(hObject, eventdata, handles)
% hObject    handle to FileVarNamePopUp (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function HausedorffZThresholdEdit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to HausedorffZThresholdEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function CurvaturePopUp_Callback(hObject, eventdata, handles)
function CurvaturePopUp_CreateFcn(hObject, eventdata, handles)
% hObject    handle to CurvaturePopUp (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
