function varargout = OFAMM_v1(varargin)
% OFAMM_v1 MATLAB code for OFAMM_v1.fig
%      OFAMM_v1, by itself, creates a new OFAMM_v1 or raises the existing
%      singleton*.
%
%      H = OFAMM_v1 returns the handle to a new OFAMM_v1 or the handle to
%      the existing singleton*.
%
%      OFAMM_v1('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in OFAMM_v1.M with the given input arguments.
%
%      OFAMM_v1('Property','Value',...) creates a new OFAMM_v1 or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before OFAMM_v1_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to OFAMM_v1_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help OFAMM_v1

% Last Modified by GUIDE v2.5 16-Jul-2019 17:56:50

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
    'gui_Singleton',  gui_Singleton, ...
    'gui_OpeningFcn', @OFAMM_v1_OpeningFcn, ...
    'gui_OutputFcn',  @OFAMM_v1_OutputFcn, ...
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
% End initialization code - DO NOT EDIT

% --- Executes just before OFAMM_v1 is made visible.
function OFAMM_v1_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to OFAMM_v1 (see VARARGIN)

% Choose default command line output for OFAMM_v1
warning('off','all')
handles.output = hObject;

if ~isempty(varargin)
    if isnumeric(varargin{1})
        handles.ImgSeq = varargin{1};
        handles.ImgSeqLoaded = 1;
        
    elseif exist(varargin{1},'dir')
        handles.PathName = varargin{1};
        FilterSpec = {'*.tif; *.raw; *.mat'};
        DialogTitle = 'Select the image sequence';
        [FileName,PathName,FilterIndex] = uigetfile(FilterSpec,DialogTitle,handles.PathName);
        if FilterIndex
            [handles.PathName, handles.FileName, handles.ExtName] = fileparts([PathName,FileName]);
            handles = loadData(handles,eventdata);
            handles.ImgSeqLoaded = 1;
        end
        
    elseif exist(varargin{1},'file')==2
        [handles.PathName, handles.FileName, handles.ExtName] = fileparts(varargin{1});
        handles = loadData(handles,eventdata);
        handles.ImgSeqLoaded = 1;
    end
end
% Update handles structure
guidata(hObject, handles);

% UIWAIT makes OFAMM_v1 wait for user response (see UIRESUME)
% uiwait(handles.figure1);

function handles = loadData(handles,eventdata)

FullFileName = fullfile(handles.PathName, [handles.FileName, handles.ExtName]);

if strcmp(handles.ExtName, '.tif') || strcmp(handles.ExtName, '.tiff')
    %     nFrames = input('Please insert number of frames for .tif file: ');
    handles.ImgSeq = imreadalltiff(FullFileName);
    handles.ImgSeqLoaded = 1;
elseif strcmp(handles.ExtName, '.raw')
    nFrames = input('Please insert number of frames for .raw file: ');
    x = input('Please insert number of pixels in "x" dir: ');
    y = input('Please insert number of pixels in "y" dir: ');
    handles.ImgSeq = imreadallraw(FullFileName,x,y,nFrames,'*float32');
    handles.ImgSeqLoaded = 1;
elseif strcmp(handles.ExtName, '.mat')
    mFile = matfile(FullFileName);
    fn = fieldnames(mFile);
    count_3d = 0;
    for Nfn = 1:length(fn)
        if numel(size(eval(['mFile.',fn{Nfn}])))==3
            count_3d = count_3d+1;
            idx_3d(count_3d) = Nfn;
        end
    end
    if ~isempty(idx_3d)
        eval(['handles.ImgSeq = mFile.',fn{idx_3d(1)},';']);
        handles.ImgSeqLoaded = 1;
    end
    delete(mFile);
else
    disp('Invalid data format. Data format could be ".tif", ".raw", or ".mat".')
end
% Load uvResults
if handles.ImgSeqLoaded == 1
    handles.uvCLGcalculated = 0;
    handles.uvHScalculated = 0;
    handles.uvTScalculated = 0;
    if exist([handles.PathName,'\uvResults.mat'],'file')
        load([handles.PathName,'\uvResults.mat']);
        if exist('uvCLG','var')
            handles.uvCLG = uvCLG;
            handles.uvCLGcalculated = 1;
            handles.FstartOFcalculated = 1;
            handles.FendOFcalculated = size(uvCLG,3);
            clear uvCLG;
        end
        if exist('uvHS','var')
            handles.uvHS = uvHS;
            handles.uvHScalculated = 1;
            handles.FstartOFcalculated = 1;
            handles.FendOFcalculated = size(uvHS,3);
            clear uvHS;
        end
        if exist('uvTS','var')
            handles.uvTS = uvTS;
            handles.uvTScalculated = 1;
            clear uvTS;
        end
        if exist('FstartOFcalculated','var')
            handles.FstartOFcalculated = FstartOFcalculated;
            clear FstartOFcalculated;
        end
        if exist('FendOFcalculated','var')
            handles.FendOFcalculated = FendOFcalculated;
            clear FendOFcalculated;
        end
    end
end
% load SourceSinkSpiralResults
if handles.ImgSeqLoaded == 1
    if exist([handles.PathName,'\SourceSinkSpiralResults.mat'],'file')
        load([handles.PathName,'\SourceSinkSpiralResults.mat']);
        if exist('SoSi','var')
            if isfield(SoSi,'CLG')
                handles.SourceCLG = SoSi.CLG.source;
                handles.SinkCLG = SoSi.CLG.sink;
                handles.ContourSourceCLG = SoSi.CLG.contour_source;
                handles.ContourSinkCLG = SoSi.CLG.contour_sink;
                handles.SSCLGcalculated=1;
            end
            if isfield(SoSi,'HS')
                handles.SourceHS = SoSi.HS.source;
                handles.SinkHS = SoSi.HS.sink;
                handles.ContourSourceHS = SoSi.HS.contour_source;
                handles.ContourSinkHS = SoSi.HS.contour_sink;
                handles.SSHScalculated=1;
            end
            clear SoSi
        end
    end
end

[dim1,dim2,nFrames] = size(handles.ImgSeq);
handles.dim1 = dim1;
handles.dim2 = dim2;
handles.nFrames = nFrames;

if isfield(handles,'Mask')
    if size(handles.Mask,1) ~= handles.dim1 || size(handles.Mask,2) ~= handles.dim2
        handles.Mask = ones(handles.dim1,handles.dim2);
        disp('The Mask size does not agrees with image sequence.\')
        disp('Mask is set to all ones.\n')
    end
else
    handles.Mask = ones(handles.dim1,handles.dim2);
end
% Get indicies inside the Mask.
[handles.rMask,handles.cMask] = find(handles.Mask > 0);
handles.idxMask = sub2ind(size(handles.Mask),handles.rMask,handles.cMask);

set(handles.FramesSlider,'Max',nFrames);
set(handles.FramesSlider,'Min',1);
minStep = 1/nFrames;
maxStep = 10/nFrames;
set(handles.FramesSlider,'SliderStep',[minStep maxStep]);
set(handles.FramesSlider,'value',1);
metaData.playFlag = 0;
set(handles.Play,'userData',metaData);
set(handles.Stop,'visible','off');
% addlistener(handles.FramesSlider,'ContinuousValueChange',@slider_frames_Callback1);
handles = loadFrame(handles,eventdata,1);

guidata(gcbo, handles);
n = 0;

function handles = loadFrame(handles,eventdata,frameNumber)
% get min and max disp values
MaxDispVal_default = str2num(get(handles.MinDispVal,'string'));
handles.MaxVal = Str2NumFromHandle(handles.MaxDispVal,MaxDispVal_default);
MinDispVal_default = str2num(get(handles.MaxDispVal,'string'));
handles.MinVal = Str2NumFromHandle(handles.MinDispVal,MinDispVal_default);
MinVal = handles.MinVal;
MaxVal = handles.MaxVal;
% imagesc frame and mask
if isfield(handles,'ImgSeq')
    if handles.ImgSeqLoaded
        if frameNumber > 0 && frameNumber <= handles.nFrames
            if get(handles.MaskChk,'value')
                thisframe = double(handles.ImgSeq(:,:,frameNumber)).*double(handles.Mask);
            else
                thisframe = handles.ImgSeq(:,:,frameNumber);
            end
            axes(handles.MainAxes);cla
            imagesc(thisframe,[MinVal,MaxVal]);
            colormap jet;
            hold on;
            set(handles.FrameNumber,'String',sprintf('%d of %d',frameNumber,handles.nFrames));
            if handles.nFrames>99
                set(handles.FrameNumber,'FontSize',8);
            end
            set(handles.FramesSlider,'value',frameNumber);
        end
    end
end
% show optical flow
if get(handles.VectorFieldVisChk,'value')
    if get(handles.CLGVisOF,'value'); methodOF = 'CLG'; end
    if get(handles.HSVisOF,'value'); methodOF = 'HS'; end
    
    if eval(['isfield(handles,','''uv',methodOF,'calculated'')'])
        if eval(['handles.uv',methodOF,'calculated'])
            if frameNumber > 0 && frameNumber <= handles.nFrames &&...
                    frameNumber >= handles.FstartOFcalculated &&...
                    frameNumber < handles.FendOFcalculated
                
                frameNumberOF = frameNumber;
                eval(['thisuv = handles.uv',methodOF,'(:,:,frameNumberOF);']);
                u = real(thisuv); if get(handles.MaskChk,'value'); u = u.*handles.Mask; end
                v = imag(thisuv); if get(handles.MaskChk,'value'); v = v.*handles.Mask; end
                [gridXDown,gridYDown,OpticalFlowDown_u] = mean_downsample(u,3);
                [~,~,OpticalFlowDown_v]                 = mean_downsample(v,3);
                axes(handles.MainAxes);
                quiver(gridXDown,gridYDown,OpticalFlowDown_u(:,:,1),OpticalFlowDown_v(:,:,1),2,'color','w','LineWidth',0.75);
                hold on;
            end
        end
    end
end

% display source/sink (Simple or Node source/sink)
markersize = 3;
linewidth = 1.2;
source_color = 'g';
sink_color = 'w';
contour_source_color = 'g';
contour_sink_color = 'w';

if get(handles.SSChk,'value')
    if get(handles.CLGVisSS,'value'); methodSS = 'CLG'; end
    if get(handles.HSVisSS,'value'); methodSS = 'HS'; end
    
    if eval(['isfield(handles,','''SS',methodSS,'calculated'')'])
        if eval(['handles.SS',methodSS,'calculated'])
            eval(['source = handles.Source',methodSS,';']);
            eval(['sink = handles.Sink',methodSS,';']);
            eval(['contour_source = handles.ContourSource',methodSS,';']);
            eval(['contour_sink = handles.ContourSink',methodSS,';']);
            %find sources and display
            [~, sourceidx] = find(source(3,:) == frameNumber);
            for idx = sourceidx
                axes(handles.MainAxes);
                plot(source(1,idx),source(2,idx),'o','color',source_color,'markersize',markersize,'markerfacecolor',source_color)
                hold on;
                plot(contour_source{idx}.xy(1,:),contour_source{idx}.xy(2,:),'-','color',contour_source_color,'linewidth',linewidth)
                hold on
            end
            %find sinks and display
            [~, sinkidx] = find(sink(3,:) == frameNumber);
            for idx = sinkidx
                axes(handles.MainAxes);
                plot(sink(1,idx),sink(2,idx),'o','color',sink_color,'markersize',markersize,'markerfacecolor',sink_color)
                hold on;
                plot(contour_sink{idx}.xy(1,:),contour_sink{idx}.xy(2,:),'-','color',contour_sink_color,'linewidth',linewidth)
                hold on
            end
        end
    end
end

% display Spiral source/sink
markersize = 3;
linewidth = 1.2;
source_color = 'b';
sink_color = 'k';
contour_source_color = 'b';
contour_sink_color = 'k';

if get(handles.SpiralChk,'value')
    if get(handles.CLGVisSpiral,'value'); methodSS = 'CLG'; end
    if get(handles.HSVisSpiral,'value'); methodSS = 'HS'; end
    
    if eval(['isfield(handles,','''SS',methodSS,'calculated'')'])
        if eval(['handles.SS',methodSS,'calculated'])
            eval(['source = handles.SourceSpiral',methodSS,';']);
            eval(['sink = handles.SinkSpiral',methodSS,';']);
            eval(['contour_source = handles.ContourSourceSpiral',methodSS,';']);
            eval(['contour_sink = handles.ContourSinkSpiral',methodSS,';']);
            %find Spiral sources and display
            [~, sourceidx] = find(source(3,:) == frameNumber);
            for idx = sourceidx
                axes(handles.MainAxes);
                plot(source(1,idx),source(2,idx),'o','color',source_color,'markersize',markersize,'markerfacecolor',source_color)
                hold on;
                plot(contour_source{idx}.xy(1,:),contour_source{idx}.xy(2,:),'-','color',contour_source_color,'linewidth',linewidth)
                hold on
            end
            %find Spiral sinks and display
            [~, sinkidx] = find(sink(3,:) == frameNumber);
            for idx = sinkidx
                axes(handles.MainAxes);
                plot(sink(1,idx),sink(2,idx),'o','color',sink_color,'markersize',markersize,'markerfacecolor',sink_color)
                hold on;
                plot(contour_sink{idx}.xy(1,:),contour_sink{idx}.xy(2,:),'-','color',contour_sink_color,'linewidth',linewidth)
                hold on
            end
        end
    end
end

% display trajectories
metaData = get(handles.Play,'userData');
if ~isfield(metaData,'playFlag')
    metaData.playFlag = 0;
    set(handles.Play,'userData',metaData);
end

if get(handles.Trajectory,'value') && ~metaData.playFlag
    handles.TrajectoryInfo.currVal = 1;
    if get(handles.CLGVisTraj,'value'); methodTraj = 'CLG'; end
    if get(handles.HSVisTraj,'value'); methodTraj = 'HS'; end
    
    if eval(['isfield(handles,','''uv',methodTraj,'calculated'')'])
        if eval(['handles.uv',methodTraj,'calculated'])
            if ~isfield(handles.TrajectoryInfo,'preVal')
                handles.TrajectoryInfo.preVal = 0;
            end
            
            if handles.TrajectoryInfo.preVal
                if handles.TrajectoryInfo.currentFrame == frameNumber
                    handles = FstartROI_fun(handles);
                    handles = FendROI_fun(handles);
                    FstartROI = max([handles.ROI.Fstart, handles.FstartOFcalculated]);
                    FendROI = min([handles.ROI.Fend,handles.FendOFcalculated]);
                    set(handles.FstartROI,'string',FstartROI);
                    set(handles.FendROI,'string',FendROI);
                else
                    FstartROI = handles.TrajectoryInfo.Fstart;
                    FstartROI = max([FstartROI, handles.FstartOFcalculated]);
                    FendROI = min([frameNumber,handles.FendOFcalculated]);
                end
                
                xy = handles.ROI.xy;
                sx = xy(1,:);
                sy = xy(2,:);
                startFrame = FstartROI;
                endFrame = FendROI;
                offsetFrame = 1;%handles.FstartOFcalculated;
                if startFrame < endFrame
                    output = Velocity_Profile(eval(['handles.uv',methodTraj]),startFrame,endFrame,offsetFrame,sx,sy);
                    % display
                    axes(handles.MainAxes);
                    str = streamline(output.str);
                    set(str,'Color','w'); hold on;
                    handles.TrajectoryInfo.preVal = 1;
                    handles.TrajectoryInfo.Fstart = FstartROI;
                end
            else
                if isfield(handles,'ROI')
                    if isfield(handles.ROI,'selected')
                        if ~handles.ROI.selected
                            handles = SelectROI(handles.SelectROI, eventdata, handles);
                        end
                    else
                        handles = SelectROI(handles.SelectROI, eventdata, handles);
                    end
                else
                    handles = SelectROI(handles.SelectROI, eventdata, handles);
                end
                
                if handles.ROI.selected
                    handles = FstartROI_fun(handles);
                    FstartROI = max([handles.ROI.Fstart, handles.FstartOFcalculated]);
                    FendROI = min([frameNumber,handles.FendOFcalculated]);
                    set(handles.FstartROI,'string',FstartROI);
                    
                    xy = handles.ROI.xy;
                    sx = xy(1,:);
                    sy = xy(2,:);
                    startFrame = FstartROI;
                    endFrame = FendROI;
                    offsetFrame = 1;
                    if startFrame < endFrame
                        output = Velocity_Profile(eval(['handles.uv',methodTraj]),startFrame,endFrame,offsetFrame,sx,sy);
                        % display
                        axes(handles.MainAxes);
                        str = streamline(output.str);
                        set(str,'Color','w'); hold on;
                        handles.TrajectoryInfo.preVal = 1;
                        handles.TrajectoryInfo.Fstart = FstartROI;
                    end
                end
            end
        end
    end
end

if get(handles.ShowROI,'Value')
    if isfield(handles,'ROI')
        if isfield(handles.ROI,'selected')
            if handles.ROI.selected
                [yROI, xROI] = find(handles.ROI.BW);
                axes(handles.MainAxes);
                plot(xROI, yROI,'m.')
            end
        end
    end
end

% --- Outputs from this function are returned to the command line.
function varargout = OFAMM_v1_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on slider movement.
function FramesSlider_Callback(hObject, eventdata, handles)
frameNumber = round(get(hObject,'Value'));
set(hObject,'Value',frameNumber);
handles = loadFrame(handles,eventdata,frameNumber);
guidata(gcbo,handles);

% --- Executes on slider movement.
function slider_frames_Callback1(hObject, eventdata, handles)

handles = guidata(gcbo);
frameNumber = round(get(handles.FramesSlider,'Value'));
set(handles.FramesSlider,'Value',frameNumber);

handles = MinDispVal_Callback(handles.MinDispVal, eventdata, handles);
handles = MaxDispVal_Callback(handles.MaxDispVal, eventdata, handles);
MinVal = handles.MinVal;
MaxVal = handles.MaxVal;

% imagesc frame and mask
if isfield(handles,'ImgSeq')
    if handles.ImgSeqLoaded
        if frameNumber > 0 && frameNumber <= handles.nFrames
            if get(handles.MaskChk,'value')
                thisframe = double(handles.ImgSeq(:,:,frameNumber)).*double(handles.Mask);
            else
                thisframe = handles.ImgSeq(:,:,frameNumber);
            end
            axes(handles.MainAxes);cla
            imagesc(thisframe,[MinVal,MaxVal]);
            colormap jet;
            hold on;
            set(handles.FrameNumber,'String',sprintf('%d of %d',frameNumber,handles.nFrames));
            set(handles.FramesSlider,'value',frameNumber);
        end
    end
end
guidata(gcbo,handles);


% --- Executes during object creation, after setting all properties.
function FramesSlider_CreateFcn(hObject, eventdata, handles)
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on button press in Play.
function Play_Callback(hObject, eventdata, handles)
set(hObject,'visible','off');
set(handles.Stop,'visible','on');
metaData = get(handles.MainAxes,'userData');
metaData.playFlag = 1;
set(hObject,'userData',metaData);
currentFrame = round(get(handles.FramesSlider,'Value'));
maxFrame = get(handles.FramesSlider,'Max');
if currentFrame == maxFrame
    set(handles.FramesSlider,'Value',1);
    currentFrame = 1;
end
while 1
    metaData = get(hObject,'userData');
    if metaData.playFlag == 0
        break;
    end
    if currentFrame <maxFrame
        currentFrame = currentFrame + 1;
        handles = loadFrame(handles,eventdata,currentFrame);
    else
        Stop_Callback(handles.Stop, eventdata, handles)
        break;
    end
    pause(0.01);
end
guidata(gcbo,handles);

% --- Executes on button press in Stop.
function Stop_Callback(hObject, eventdata, handles)
metaData = get(handles.Play,'userData');
metaData.playFlag = 0;
set(handles.Play,'userData',metaData);
set(hObject,'visible','off');
set(handles.Play,'visible','on');
pause(0.3);
guidata(gcbo,handles);

% --- Executes when selected object is changed in OFParameters.
function OFParameters_SelectionChangeFcn(hObject, eventdata, handles)
currentFrame = round(get(handles.FramesSlider,'Value'));
handles = loadFrame(handles,eventdata,currentFrame);
guidata(gcbo,handles);

% --- Executes on button press in VectorFieldVisChk.
function VectorFieldVisChk_Callback(hObject, eventdata, handles)
currentFrame = round(get(handles.FramesSlider,'Value'));
handles = loadFrame(handles,eventdata,currentFrame);
guidata(gcbo,handles);

% --- Executes on button press in pushbutton_selectPoints.
function pushbutton_selectPoints_Callback(hObject, eventdata, handles)
guidata(gcbo,handles);

function edit_pointsToTrack_Callback(hObject, eventdata, handles)
guidata(gcbo,handles);

% --- Executes during object creation, after setting all properties.
function edit_pointsToTrack_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function alphaCLG_Callback(hObject, eventdata, handles)

% --- Executes during object creation, after setting all properties.
function alphaCLG_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function ratioCLG_Callback(hObject, eventdata, handles)

% --- Executes during object creation, after setting all properties.
function ratioCLG_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function minWidthCLG_Callback(hObject, eventdata, handles)

% --- Executes during object creation, after setting all properties.
function minWidthCLG_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function nOuterFPIterationsCLG_Callback(hObject, eventdata, handles)


% --- Executes during object creation, after setting all properties.
function nOuterFPIterationsCLG_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function nInnerFPIterationsCLG_Callback(hObject, eventdata, handles)

% --- Executes during object creation, after setting all properties.
function nInnerFPIterationsCLG_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function nSORIterationsCLG_Callback(hObject, eventdata, handles)

% --- Executes during object creation, after setting all properties.
function nSORIterationsCLG_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function alphaHS_Callback(hObject, eventdata, handles)

% --- Executes during object creation, after setting all properties.
function alphaHS_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function IterationsHS_Callback(hObject, eventdata, handles)

% --- Executes during object creation, after setting all properties.
function IterationsHS_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function CorrWinPixelsTS_Callback(hObject, eventdata, handles)

% --- Executes during object creation, after setting all properties.
function CorrWinPixelsTS_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function CorrWinFramesTS_Callback(hObject, eventdata, handles)

% --- Executes during object creation, after setting all properties.
function CorrWinFramesTS_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function TargetFramesTS_Callback(hObject, eventdata, handles)

% --- Executes during object creation, after setting all properties.
function TargetFramesTS_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in RunOF.
function RunOF_Callback(hObject, eventdata, handles)
handles = RunOpticalFlowAnalysisButton(handles);
guidata(gcbo,handles);
n=0;


% --- Executes on button press in RunSS.
function RunSS_Callback(hObject, eventdata, handles)
if isfield(handles,'uvCLGcalculated')
    if handles.uvCLGcalculated
        method = 'CLG';
        handles = calc_save_SourceSink(handles,method);
    end
end
if isfield(handles,'uvHScalculated')
    if handles.uvHScalculated
        method = 'HS';
        handles = calc_save_SourceSink(handles,method);
    end
end
guidata(gcbo,handles);

% --- Executes on button press in ROISpeedAngle.
function ROISpeedAngle_Callback(hObject, eventdata, handles)

% --- Executes on button press in ViewInstSpeedAngle.
function ViewInstSpeedAngle_Callback(hObject, eventdata, handles)
handles = InstSpDir(handles);
if isfield(handles,'ROI')
    if isfield(handles.ROI,'selected')
        if handles.ROI.selected
            plotflag = 1;
            handles = plotInstSpDir(handles, plotflag);
            handles.InstSpAngcalculated = 1;
        else
            handles.InstSpAngcalculated = 0;
        end
    else
        handles.InstSpAngcalculated = 0;
    end
else
    handles.InstSpAngcalculated = 0;
end

guidata(gcbo,handles);


% --- Executes on button press in pushbutton11.
function pushbutton11_Callback(hObject, eventdata, handles)

% --- Executes on button press in checkbox5.
function checkbox5_Callback(hObject, eventdata, handles)

% --- Executes on button press in LoadImgSeq.
function LoadImgSeq_Callback(hObject, eventdata, handles)
FilterSpec = {'*.tif; *.raw; *.mat'};
DialogTitle = 'Select the image sequence';
[FileName,PathName,FilterIndex] = uigetfile(FilterSpec,DialogTitle);
if FilterIndex
    [handles.PathName, handles.FileName, handles.ExtName] = fileparts([PathName,FileName]);
    handles = loadData(handles,eventdata);
    handles.ImgSeqLoaded = 1;
end
guidata(hObject, handles);

% --- Executes on button press in LoadMask.
function LoadMask_Callback(hObject, eventdata, handles)
FilterSpec = {'*.tif; *.raw; *.mat'};
DialogTitle = 'Select the mask image';
[MaskFileName,MaskPathName,FilterIndex] = uigetfile(FilterSpec,DialogTitle);
if FilterIndex
    [handles.MaskPathName, handles.MaskFileName, handles.MaskExtName] = fileparts([MaskPathName,MaskFileName]);
    handles = loadMaskData(handles);
    handles.Mask = handles.Mask > 0; % Make sure the maske is binary
    % Get indicies inside the Mask.
    [handles.rMask,handles.cMask] = find(handles.Mask > 0);
    handles.idxMask = sub2ind(size(handles.Mask),handles.rMask,handles.cMask);
    guidata(hObject, handles);
end
n = 0;

function handles = loadMaskData(handles)
FullMaskFileName = fullfile(handles.MaskPathName, [handles.MaskFileName, handles.MaskExtName]);

if strcmp(handles.MaskExtName, '.tif') || strcmp(handles.MaskExtName, '.tiff')
    handles.Mask = imread(FullMaskFileName);
elseif strcmp(handles.MaskExtName, '.raw')
    x = input('Please insert number of pixels in "x" dir: ');
    y = input('Please insert number of pixels in "y" dir: ');
    handles.Mask = imreadallraw(FullMaskFileName,x,y,1,'*float32');
elseif strcmp(handles.MaskExtName, '.mat')
    mFile = matfile(FullMaskFileName);
    fn = fieldnames(mFile);
    Nfn = 1;
    keepgoing = 1;
    while keepgoing
        a = eval(['mFile.',fn{Nfn}]);
        if isnumeric(a)
            if size(a,3) == 1
                eval(['handles.Mask = mFile.',fn{Nfn},';']);
                handles.ImgSeqLoaded = 1;
                keepgoing = 0;
            end
        end
        Nfn = Nfn+1;
        if Nfn > length(fn)
            keepgoing = 0;
        end
    end
    delete(mFile);
else
    handles.Mask = ones(handles.dim1,handles.dim2);
    disp('Invalid data format. Data format could be ".tif", ".raw", or ".mat".\n')
    disp('Mask is set to all ones.\n')
end
if isfield(handles,'dim1') && isfield(handles,'dim2')
    if size(handles.Mask,1) ~= handles.dim1 || size(handles.Mask,2) ~= handles.dim2
        handles.Mask = ones(handles.dim1,handles.dim2);
        disp('The Mask size does not agrees with image sequence.\')
        disp('Mask is set to all ones.\n')
    end
end


function MaxLagTS_Callback(hObject, eventdata, handles)

% --- Executes during object creation, after setting all properties.
function MaxLagTS_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in MaskChk.
function MaskChk_Callback(hObject, eventdata, handles)
currentFrame = round(get(handles.FramesSlider,'Value'));
handles = loadFrame(handles,eventdata,currentFrame);
guidata(gcbo,handles);
n=0;

% --- Executes on button press in IntensityTime.
function handles = IntensityTime_Callback(hObject, eventdata, handles)
[~, flag] = getimage(handles.MainAxes); % check if there is an image in the main axes
if flag
    if isfield(handles,'ROI')
        if isfield(handles.ROI,'selected')
            if ~handles.ROI.selected
                handles = SelectROI_Callback(handles.SelectROI, eventdata, handles);
            end
        else
            handles = SelectROI_Callback(handles.SelectROI, eventdata, handles);
        end
    else
        handles = SelectROI_Callback(handles.SelectROI, eventdata, handles);
    end
    
    handles = FstartROI_fun(handles);
    handles = FendROI_fun(handles);
    FstartROI = handles.ROI.Fstart;
    FendROI = handles.ROI.Fend;
    handles.IntensitySig = zeros(1,FendROI-FstartROI+1);
    for idx = FstartROI:FendROI
        thisframe = handles.ImgSeq(:,:,idx);
        handles.IntensitySig(idx-FstartROI+1) = mean(thisframe(handles.ROI.BW));
    end
    if isfield(handles,'IntensityFig')
        if ishandle(handles.IntensityFig)
            figure(handles.IntensityFig);
            cla;plot(FstartROI:FendROI, handles.IntensitySig)
            xlabel('frame number')
            ylabel('Mean Intensity')
        else
            %             figure(101,'visible','off');
            figure(101);
            handles.IntensityFig = gcf;
            set(handles.IntensityFig,'visible','on','numbertitle','off','name','Intensity vs. Time for selected ROI');
            plot(FstartROI:FendROI, handles.IntensitySig)
            xlabel('frame number')
            ylabel('Mean Intensity')
        end
    else
        figure(101);
        handles.IntensityFig = gcf;
        set(handles.IntensityFig,'visible','on','numbertitle','off','name','Intensity vs. Time for selected ROI');
        plot(FstartROI:FendROI,handles.IntensitySig)
        xlabel('frame number')
        ylabel('Mean Intensity')
    end
end
guidata(gcbo,handles);
n=0;


% --- Executes on button press in ViewTempSpeedTrajLength.
function ViewTempSpeedTrajLength_Callback(hObject, eventdata, handles)
handles = TempSpTrajLen(handles);
if isfield(handles,'ROI')
    if isfield(handles.ROI,'selected')
        if handles.ROI.selected
            plotflag = 1;
            handles = plotTempSpTrajLen(handles, plotflag);
            handles.TempSpTrajLengthcalculated = 1;
        else
            handles.TempSpTrajLengthcalculated = 0;
        end
    else
        handles.TempSpTrajLengthcalculated = 0;
    end
else
    handles.TempSpTrajLengthcalculated = 0;
end
guidata(gcbo,handles);
n=0;

% --- Executes on button press in SaveSpeedAngle.
function SaveSpeedAngle_Callback(hObject, eventdata, handles)
handles = SaveSpeedAngle_fun(handles);
guidata(gcbo,handles);


% --- Executes on button press in ViewTraj.
function ViewTraj_Callback(hObject, eventdata, handles)


% --- Executes during object creation, after setting all properties.
function ROISaveNameSpeedAngle_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in SSChk.
function SSChk_Callback(hObject, eventdata, handles)
currentFrame = round(get(handles.FramesSlider,'Value'));
handles = loadFrame(handles,eventdata,currentFrame);
guidata(gcbo,handles);

% --- Executes on button press in SaveSS.
function SaveSS_Callback(hObject, eventdata, handles)


% --- Executes on button press in SaveOF.
function SaveOF_Callback(hObject, eventdata, handles)


% --- Executes on button press in HSsave.
function HSsave_Callback(hObject, eventdata, handles)


% --- Executes on button press in CLGsave.
function CLGsave_Callback(hObject, eventdata, handles)


% --- Executes on button press in TSsave.
function TSsave_Callback(hObject, eventdata, handles)


% --- Executes on button press in OFrawFile.
function OFrawFile_Callback(hObject, eventdata, handles)


% --- Executes on button press in OFmatFile.
function OFmatFile_Callback(hObject, eventdata, handles)



function FstartOF_Callback(hObject, eventdata, handles)


% --- Executes during object creation, after setting all properties.
function FstartOF_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function FendOF_Callback(hObject, eventdata, handles)


% --- Executes during object creation, after setting all properties.
function FendOF_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function handles = MinDispVal_Callback(hObject, eventdata, handles)
MinDispVal_default = str2num(get(handles.MaxDispVal,'string'));
handles.MinVal = Str2NumFromHandle(hObject,MinDispVal_default);
currentFrame = round(get(handles.FramesSlider,'Value'));
handles = loadFrame(handles,eventdata,currentFrame);
guidata(gcbo,handles);


% --- Executes during object creation, after setting all properties.
function MinDispVal_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function handles = MaxDispVal_Callback(hObject, eventdata, handles)
MaxDispVal_default = str2num(get(handles.MinDispVal,'string'));
handles.MaxVal = Str2NumFromHandle(hObject,MaxDispVal_default);
currentFrame = round(get(handles.FramesSlider,'Value'));
handles = loadFrame(handles,eventdata,currentFrame);
guidata(gcbo,handles);


% --- Executes during object creation, after setting all properties.
function MaxDispVal_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in SelectROI.
function handles = SelectROI_Callback(hObject, eventdata, handles)
handles = SelectROI(hObject, eventdata, handles);
guidata(gcbo,handles);


% --- Executes on button press in SelectPoints.
function handles = SelectPoints_Callback(hObject, eventdata, handles)
[~, flag] = getimage(handles.MainAxes); % check if there is an image in the main axes
if flag
    [x, y] = ginput;
    x = round(x);
    y = round(y);
    InRange = (x >= 1) & (x <= handles.dim2) & (y >= 1) & (y <= handles.dim1);
    x = x(InRange);
    y = y(InRange);
    handles.ROI.xy = [x';y'];
    handles.ROI.BW = zeros(handles.dim1,handles.dim2)>0;
    for p = 1:length(x)
        handles.ROI.BW(y(p),x(p)) = 1;
    end
    handles.ROI.selected = 1;
    guidata(gcbo,handles);
end


% --- Executes during object creation, after setting all properties.
function FstartROI_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes during object creation, after setting all properties.
function FendROI_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in checkbox16.
function checkbox16_Callback(hObject, eventdata, handles)


% --- Executes on button press in Trajectory.
function Trajectory_Callback(hObject, eventdata, handles)

if get(hObject,'value')
    handles.TrajectoryInfo.currVal = 1;
    if get(handles.CLGVisTraj,'value'); methodTraj = 'CLG'; end
    if get(handles.HSVisTraj,'value'); methodTraj = 'HS'; end
    
    handles.TrajectoryInfo.currentFrame = round(get(handles.FramesSlider,'Value'));
    
    if eval(['isfield(handles,','''uv',methodTraj,'calculated'')'])
        if eval(['handles.uv',methodTraj,'calculated'])
            
            if isfield(handles,'ROI')
                if isfield(handles.ROI,'selected')
                    if ~handles.ROI.selected
                        handles = SelectROI(handles.SelectROI, eventdata, handles);
                    end
                else
                    handles = SelectROI(handles.SelectROI, eventdata, handles);
                end
            else
                handles = SelectROI(handles.SelectROI, eventdata, handles);
            end
            
            if handles.ROI.selected
                handles = FstartROI_fun(handles);
                handles = FendROI_fun(handles);
                FstartROI = max([handles.ROI.Fstart, handles.FstartOFcalculated]);
                FendROI = min([handles.ROI.Fend,handles.FendOFcalculated]);
                set(handles.FstartROI,'string',FstartROI);
                set(handles.FendROI,'string',FendROI);
                
                xy = handles.ROI.xy;
                sx = xy(1,:);
                sy = xy(2,:);
                startFrame = FstartROI;
                endFrame = FendROI;
                offsetFrame = 1;
                if startFrame < endFrame
                    output = Velocity_Profile(eval(['handles.uv',methodTraj]),startFrame,endFrame,offsetFrame,sx,sy);
                    % display
                    str = streamline(output.str);
                    set(str,'Color','w'); hold on;
                    handles.TrajectoryInfo.preVal = 1;
                    handles.TrajectoryInfo.Fstart = FstartROI;
                end
            end
        end
    end
else
    currentFrame = round(get(handles.FramesSlider,'Value'));
    handles = loadFrame(handles,eventdata,currentFrame);
end

if ~get(hObject,'value')
    handles.TrajectoryInfo.preVal = 0;
    handles.TrajectoryInfo.Fstart = 0;
end
guidata(gcbo,handles);


% --- Executes during object creation, after setting all properties.
function NbinSpDir_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function handles = NbinSpDir_Callback(hObject, eventdata, handles)
handles = NbinSpDir_fun(handles);
guidata(gcbo,handles);


function handles = NbinSpTraj_Callback(hObject, eventdata, handles)
handles = NbinSpTraj_fun(handles);
guidata(gcbo,handles);


% --- Executes during object creation, after setting all properties.
function NbinSpTraj_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function handles = ROISaveNameSpeedAngle_Callback(hObject, eventdata, handles)
handles = ROISaveNameSpeedAngle_fun(handles);
guidata(gcbo,handles);

function handles = FendROI_Callback(hObject, eventdata, handles)
handles = FendROI_fun(handles);
guidata(gcbo,handles);


function handles = FstartROI_Callback(hObject, eventdata, handles)
handles = FstartROI_fun(handles);
guidata(gcbo,handles);


% --- Executes on button press in CLGVisSS.
function CLGVisSS_Callback(hObject, eventdata, handles)
if get(hObject,'value')
    set(handles.HSVisSS,'value',0);
else
    set(handles.HSVisSS,'value',1);
end
currentFrame = round(get(handles.FramesSlider,'Value'));
handles = loadFrame(handles,eventdata,currentFrame);
guidata(gcbo,handles);

% --- Executes on button press in HSVisSS.
function HSVisSS_Callback(hObject, eventdata, handles)
if get(hObject,'value')
    set(handles.CLGVisSS,'value',0);
else
    set(handles.CLGVisSS,'value',1);
end
currentFrame = round(get(handles.FramesSlider,'Value'));
handles = loadFrame(handles,eventdata,currentFrame);
guidata(gcbo,handles);

% --- Executes on button press in CLGVisTraj.
function CLGVisTraj_Callback(hObject, eventdata, handles)
if get(hObject,'value')
    set(handles.HSVisTraj,'value',0);
else
    set(handles.HSVisTraj,'value',1);
end
% Trajectory_Callback(handles.Trajectory, eventdata, handles);
currentFrame = round(get(handles.FramesSlider,'Value'));
handles = loadFrame(handles,eventdata,currentFrame);
guidata(gcbo,handles);

% --- Executes on button press in HSVisTraj.
function HSVisTraj_Callback(hObject, eventdata, handles)
if get(hObject,'value')
    set(handles.CLGVisTraj,'value',0);
else
    set(handles.CLGVisTraj,'value',1);
end
% Trajectory_Callback(handles.Trajectory, eventdata, handles);
currentFrame = round(get(handles.FramesSlider,'Value'));
handles = loadFrame(handles,eventdata,currentFrame);
guidata(gcbo,handles);

% --- Executes on button press in CLGVisOF.
function CLGVisOF_Callback(hObject, eventdata, handles)
if get(hObject,'value')
    set(handles.HSVisOF,'value',0);
else
    set(handles.HSVisOF,'value',1);
end
VectorFieldVisChk_Callback(handles.VectorFieldVisChk, eventdata, handles);
guidata(gcbo,handles);

% --- Executes on button press in HSVisOF.
function HSVisOF_Callback(hObject, eventdata, handles)
if get(hObject,'value')
    set(handles.CLGVisOF,'value',0);
else
    set(handles.CLGVisOF,'value',1);
end
VectorFieldVisChk_Callback(handles.VectorFieldVisChk, eventdata, handles);
guidata(gcbo,handles);


% --- Executes on button press in runCLG.
function runCLG_Callback(hObject, eventdata, handles)


% --- Executes on button press in runHS.
function runHS_Callback(hObject, eventdata, handles)


% --- Executes on button press in runTS.
function runTS_Callback(hObject, eventdata, handles)


% --- Executes on button press in SpiralChk.
function SpiralChk_Callback(hObject, eventdata, handles)
currentFrame = round(get(handles.FramesSlider,'Value'));
handles = loadFrame(handles,eventdata,currentFrame);
guidata(gcbo,handles);

% --- Executes on button press in CLGVisSpiral.
function CLGVisSpiral_Callback(hObject, eventdata, handles)
if get(hObject,'value')
    set(handles.HSVisSpiral,'value',0);
else
    set(handles.HSVisSpiral,'value',1);
end
currentFrame = round(get(handles.FramesSlider,'Value'));
handles = loadFrame(handles,eventdata,currentFrame);
guidata(gcbo,handles);

% --- Executes on button press in HSVisSpiral.
function HSVisSpiral_Callback(hObject, eventdata, handles)
if get(hObject,'value')
    set(handles.CLGVisSpiral,'value',0);
else
    set(handles.CLGVisSpiral,'value',1);
end
currentFrame = round(get(handles.FramesSlider,'Value'));
handles = loadFrame(handles,eventdata,currentFrame);
guidata(gcbo,handles);


% --- Executes on button press in ZoomChk.
function ZoomChk_Callback(hObject, eventdata, handles)
axes(handles.MainAxes);
if get(hObject,'Value')
    zoom on
else
    zoom off
end
guidata(gcbo,handles);

% --- Executes on button press in CursorChk.
function CursorChk_Callback(hObject, eventdata, handles)
axes(handles.MainAxes);
if get(hObject,'Value')
    datacursormode on
else
    datacursormode off
end
guidata(gcbo,handles);

% --- Executes on button press in Help.
function Help_Callback(hObject, eventdata, handles)
OFAMMfilePath = mfilename('fullpath');
slashPos = find(OFAMMfilePath == '\');
PathName = OFAMMfilePath(1:slashPos(end));
ReadmePathName = [PathName 'Readme.txt'];
eval(['!notepad ' ReadmePathName])
guidata(gcbo,handles);

% --- Executes on button press in AboutOFAMM.
function AboutOFAMM_Callback(hObject, eventdata, handles)
if ~isfield(handles,'AboutFig')
    posMonitor = get(0,'MonitorPositions');
    w = posMonitor(3)/6;
    posAboutFig(3) = w;
    posAboutFig(4) = w;%*posMonitor(3)/posMonitor(4);
    posAboutFig(1) = posMonitor(3)/2-posAboutFig(3)/2;
    posAboutFig(2) = posMonitor(4)/2-0*posAboutFig(4)/2;
    figure('MenuBar','none','ToolBar','none','units',get(0,'units'),'position',posAboutFig);
    
    handles.AboutFig = gcf;
    set(handles.AboutFig,'visible','on','numbertitle','off','Resize','on','name','About OFAMM');
end
AboutOFAMMWin
guidata(gcbo,handles);

function AboutOFAMMWin
a = 0.2;
x = -3:a:3;
y = -3:a:3;
ax=axes;
set(ax,'position',[0 0 1 1])
axis off
[xx,yy] = meshgrid(x,y);
zz = peaks(xx,yy);
hold on
pcolor(x,y,zz);
axis([-3 3 -3 3]);
colormap((jet+white)/2);
shading interp
[px,py] = gradient(zz,.2,.2);
c = [1 1 1]*0.7;
quiver(x,y,px,py,2,'color',c);

maxX = max(x); minX = min(x); dx = maxX - minX;
maxY = max(y); minY = min(y); dy = maxY - minY;

txtStr = 'OFAMM  v.1.0';
xt = minX + dx*0.5;
yt = minY + dy*0.9;
text(xt,yt,txtStr,'fontsize',14,'FontWeight','bold','fontname','times','HorizontalAlignment','center');

txtStr = 'A toolbox to investigate the spatiotemporal';
xt = minX + dx*0.5;
yt = minY + dy*0.8;
text(xt,yt,txtStr,'fontsize',8,'FontWeight','bold','fontname','arial','HorizontalAlignment','center');

txtStr = 'dynamics of mesoscale brain activity.';
xt = minX + dx*0.5;
yt = minY + dy*0.75;
text(xt,yt,txtStr,'fontsize',8,'FontWeight','bold','fontname','arial','HorizontalAlignment','center');

txtStr = 'by:';
xt = minX + dx*0.5;
yt = minY + dy*0.65;
text(xt,yt,txtStr,'fontsize',9,'FontWeight','bold','fontname','times','HorizontalAlignment','center');

a = 0.6; b = 0.08;
txtStr = 'Navvab Afrashteh';
xt = minX + dx*0.5;
yt = minY + dy*(a-0*b);
text(xt,yt,txtStr,'fontsize',9,'FontWeight','bold','fontname','times','HorizontalAlignment','center');

txtStr = 'Samsoon Inayat';
xt = minX + dx*0.5;
yt = minY + dy*(a-1*b);
text(xt,yt,txtStr,'fontsize',9,'FontWeight','bold','fontname','times','HorizontalAlignment','center');

txtStr = 'Mostafa Mohsenvand';
xt = minX + dx*0.5;
yt = minY + dy*(a-2*b);
text(xt,yt,txtStr,'fontsize',9,'FontWeight','bold','fontname','times','HorizontalAlignment','center');

txtStr = 'Majid H. Mohajerani';
xt = minX + dx*0.5;
yt = minY + dy*(a-3*b);
text(xt,yt,txtStr,'fontsize',9,'FontWeight','bold','fontname','times','HorizontalAlignment','center');

txtStr = 'Canadian Centre for Behavioural Neuroscience';
xt = minX + dx*0.5;
yt = minY + dy*(a-4.25*b);
text(xt,yt,txtStr,'fontsize',9,'FontWeight','bold','fontname','times','HorizontalAlignment','center');

txtStr = 'University of Lethbridge';
xt = minX + dx*0.5;
yt = minY + dy*(a-5*b);
text(xt,yt,txtStr,'fontsize',9,'FontWeight','bold','fontname','times','HorizontalAlignment','center');

txtStr = 'For more information visit:';
xt = minX + dx*0.5;
yt = minY + dy*(a-6*b);
text(xt,yt,txtStr,'fontsize',10,'FontWeight','normal','fontname','times','HorizontalAlignment','center');

cbStr = 'web(''http://lethbridgebraindynamics.com/OFAMM/'');';
txtStr = 'Lethbridge Brain Dynamics';
xt = minX + dx*0.5;
yt = minY + dy*(a-6.7*b);
htxt = text(xt,yt,txtStr,'fontsize',10,'FontWeight','normal','fontname','times','HorizontalAlignment','center','color','b');
set(htxt,'ButtonDownFcn',cbStr)
hold off


% --- Executes on button press in ShowROI.
function ShowROI_Callback(hObject, eventdata, handles)
currentFrame = round(get(handles.FramesSlider,'Value'));
handles = loadFrame(handles,eventdata,currentFrame);
guidata(gcbo,handles);
