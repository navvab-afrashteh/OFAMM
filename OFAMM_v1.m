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

% Last Modified by GUIDE v2.5 09-Feb-2017 18:42:29

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
% handles = MinDispVal_Callback(handles.MinDispVal, eventdata, handles);
% handles = MaxDispVal_Callback(handles.MaxDispVal, eventdata, handles);
MaxDispVal_default = 1;
handles.MaxVal = Str2NumFromHandle(handles.MaxDispVal,MaxDispVal_default);
MinDispVal_default = 0;
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
                
                frameNumberOF = frameNumber-handles.FstartOFcalculated+1;
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
%                     FstartROI = handles.TrajectoryInfo.Fstart;
%                     FstartROI = max([FstartROI, handles.FstartOFcalculated]);
%                     FendROI = min([frameNumber,handles.FendOFcalculated]);
                    
                    handles = FstartROI_Callback(handles.FstartROI, eventdata, handles);
                    handles = FendROI_Callback(handles.FendROI, eventdata, handles);
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
                startFrame = FstartROI - handles.FstartOFcalculated +1;
                endFrame = FendROI - handles.FstartOFcalculated +1;
                offsetFrame = handles.FstartOFcalculated;
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
                    handles = FstartROI_Callback(handles.FstartROI, eventdata, handles);
                    FstartROI = max([handles.ROI.Fstart, handles.FstartOFcalculated]);
                    FendROI = min([frameNumber,handles.FendOFcalculated]);
                    set(handles.FstartROI,'string',FstartROI);
                    
                    xy = handles.ROI.xy;
                    sx = xy(1,:);
                    sy = xy(2,:);
                    startFrame = FstartROI - handles.FstartOFcalculated +1;
                    endFrame = FendROI - handles.FstartOFcalculated +1;
                    offsetFrame = handles.FstartOFcalculated;
                    if startFrame < endFrame
                        output = Velocity_Profile(eval(['handles.uv',methodTraj]),startFrame,endFrame,offsetFrame,sx,sy);
                        % display
                        axes(handles.MainAxes);
                        str = streamline(output.str);
                        set(str,'Color','r');
                        hold on
                        handles.TrajectoryInfo.preVal = 1;
                        handles.TrajectoryInfo.Fstart = FstartROI;
                    end
                end
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
% hObject    handle to FramesSlider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
frameNumber = round(get(hObject,'Value'));
set(hObject,'Value',frameNumber);
handles = loadFrame(handles,eventdata,frameNumber);
guidata(gcbo,handles);

% --- Executes on slider movement.
function slider_frames_Callback1(hObject, eventdata, handles)
% hObject    handle to FramesSlider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
% handles = guihandles;
% hObject = handles.FramesSlider;
% frameNumber = round(get(hObject,'Value'));
% set(hObject,'Value',frameNumber);
% handles = guidata(gcbo);
% if frameNumber > 0
%     thisFrame = handles.ImgSeq(:,:,frameNumber);
% end
% imagesc(thisFrame,'Parent',handles.MainAxes);
% colormap jet;
% set(handles.FrameNumber,'String',sprintf('%d of %d',frameNumber,handles.nFrames));
% set(handles.FramesSlider,'value',frameNumber);
% set(handles.MainAxes,'userdata',frameNumber);

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
%             set(handles.MainAxes,'userdata',frameNumber);
        end
    end
end


% --- Executes during object creation, after setting all properties.
function FramesSlider_CreateFcn(hObject, eventdata, handles)
% hObject    handle to FramesSlider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on button press in Play.
function Play_Callback(hObject, eventdata, handles)
% hObject    handle to Play (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
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
%         pushbutton_stop_Callback(handles.Stop, eventdata, handles)
        Stop_Callback(handles.Stop, eventdata, handles)
        break;
    end
    pause(0.01);
end

% --- Executes on button press in Stop.
function Stop_Callback(hObject, eventdata, handles)
% hObject    handle to Stop (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
metaData = get(handles.Play,'userData');
metaData.playFlag = 0;
set(handles.Play,'userData',metaData);
set(hObject,'visible','off');
set(handles.Play,'visible','on');
pause(0.3);


% --- Executes on button press in pushbutton_selectROI.
% function pushbutton_selectROI_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_selectROI (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% BW = roipoly;
% [rows cols] = find(BW);
% handles.mask = BW;
% checkbox_trackROI_Callback(handles.checkbox_trackROI, eventdata, handles);
% guidata(gcbo,handles);

% --- Executes on button press in checkbox_trackROI.
% function checkbox_trackROI_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_trackROI (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_trackROI

% frameNumber = round(get(handles.FramesSlider,'Value'));
% set(hObject,'userdata',frameNumber);

% --- Executes when selected object is changed in OFParameters.
function OFParameters_SelectionChangeFcn(hObject, eventdata, handles)
% hObject    handle to the selected object in OFParameters 
% eventdata  structure with the following fields (see UIBUTTONGROUP)
%	EventName: string 'SelectionChanged' (read only)
%	OldValue: handle of the previously selected object or empty if none was selected
%	NewValue: handle of the currently selected object
% handles    structure with handles and user data (see GUIDATA)
currentFrame = round(get(handles.FramesSlider,'Value'));
handles = loadFrame(handles,eventdata,currentFrame);

% --- Executes on button press in VectorFieldVisChk.
function VectorFieldVisChk_Callback(hObject, eventdata, handles)
% hObject    handle to VectorFieldVisChk (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of VectorFieldVisChk
currentFrame = round(get(handles.FramesSlider,'Value'));
handles = loadFrame(handles,eventdata,currentFrame);
guidata(gcbo,handles);

% --- Executes on button press in pushbutton_selectPoints.
function pushbutton_selectPoints_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_selectPoints (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% [x y] = ginput;
% col = round(x);
% row = round(y);
% mask = zeros(size(handles.ImgSeq(:,:,1)));
% for ii = 1:length(row)
%        mask(row(ii),col(ii)) = 1;
% end
% handles.mask = mask;
% checkbox_trackROI_Callback(handles.checkbox_trackROI, eventdata, handles);
% guidata(gcbo,handles);
% set(handles.edit_pointsToTrack,'String',sprintf('%d,%d',row,col));
% nothin = 0;



function edit_pointsToTrack_Callback(hObject, eventdata, handles)
% hObject    handle to edit_pointsToTrack (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_pointsToTrack as text
%        str2double(get(hObject,'String')) returns contents of edit_pointsToTrack as a double


% --- Executes during object creation, after setting all properties.
function edit_pointsToTrack_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_pointsToTrack (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function alphaCLG_Callback(hObject, eventdata, handles)
% hObject    handle to alphaCLG (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of alphaCLG as text
%        str2double(get(hObject,'String')) returns contents of alphaCLG as a double


% --- Executes during object creation, after setting all properties.
function alphaCLG_CreateFcn(hObject, eventdata, handles)
% hObject    handle to alphaCLG (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function ratioCLG_Callback(hObject, eventdata, handles)
% hObject    handle to ratioCLG (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of ratioCLG as text
%        str2double(get(hObject,'String')) returns contents of ratioCLG as a double


% --- Executes during object creation, after setting all properties.
function ratioCLG_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ratioCLG (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function minWidthCLG_Callback(hObject, eventdata, handles)
% hObject    handle to minWidthCLG (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of minWidthCLG as text
%        str2double(get(hObject,'String')) returns contents of minWidthCLG as a double


% --- Executes during object creation, after setting all properties.
function minWidthCLG_CreateFcn(hObject, eventdata, handles)
% hObject    handle to minWidthCLG (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function nOuterFPIterationsCLG_Callback(hObject, eventdata, handles)
% hObject    handle to nOuterFPIterationsCLG (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of nOuterFPIterationsCLG as text
%        str2double(get(hObject,'String')) returns contents of nOuterFPIterationsCLG as a double


% --- Executes during object creation, after setting all properties.
function nOuterFPIterationsCLG_CreateFcn(hObject, eventdata, handles)
% hObject    handle to nOuterFPIterationsCLG (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function nInnerFPIterationsCLG_Callback(hObject, eventdata, handles)
% hObject    handle to nInnerFPIterationsCLG (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of nInnerFPIterationsCLG as text
%        str2double(get(hObject,'String')) returns contents of nInnerFPIterationsCLG as a double


% --- Executes during object creation, after setting all properties.
function nInnerFPIterationsCLG_CreateFcn(hObject, eventdata, handles)
% hObject    handle to nInnerFPIterationsCLG (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function nSORIterationsCLG_Callback(hObject, eventdata, handles)
% hObject    handle to nSORIterationsCLG (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of nSORIterationsCLG as text
%        str2double(get(hObject,'String')) returns contents of nSORIterationsCLG as a double


% --- Executes during object creation, after setting all properties.
function nSORIterationsCLG_CreateFcn(hObject, eventdata, handles)
% hObject    handle to nSORIterationsCLG (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function alphaHS_Callback(hObject, eventdata, handles)
% hObject    handle to alphaHS (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of alphaHS as text
%        str2double(get(hObject,'String')) returns contents of alphaHS as a double


% --- Executes during object creation, after setting all properties.
function alphaHS_CreateFcn(hObject, eventdata, handles)
% hObject    handle to alphaHS (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function IterationsHS_Callback(hObject, eventdata, handles)
% hObject    handle to IterationsHS (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of IterationsHS as text
%        str2double(get(hObject,'String')) returns contents of IterationsHS as a double


% --- Executes during object creation, after setting all properties.
function IterationsHS_CreateFcn(hObject, eventdata, handles)
% hObject    handle to IterationsHS (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function CorrWinPixelsTS_Callback(hObject, eventdata, handles)
% hObject    handle to CorrWinPixelsTS (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of CorrWinPixelsTS as text
%        str2double(get(hObject,'String')) returns contents of CorrWinPixelsTS as a double


% --- Executes during object creation, after setting all properties.
function CorrWinPixelsTS_CreateFcn(hObject, eventdata, handles)
% hObject    handle to CorrWinPixelsTS (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function CorrWinFramesTS_Callback(hObject, eventdata, handles)
% hObject    handle to CorrWinFramesTS (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of CorrWinFramesTS as text
%        str2double(get(hObject,'String')) returns contents of CorrWinFramesTS as a double


% --- Executes during object creation, after setting all properties.
function CorrWinFramesTS_CreateFcn(hObject, eventdata, handles)
% hObject    handle to CorrWinFramesTS (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function TargetFramesTS_Callback(hObject, eventdata, handles)
% hObject    handle to TargetFramesTS (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of TargetFramesTS as text
%        str2double(get(hObject,'String')) returns contents of TargetFramesTS as a double


% --- Executes during object creation, after setting all properties.
function TargetFramesTS_CreateFcn(hObject, eventdata, handles)
% hObject    handle to TargetFramesTS (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in RunOF.
function RunOF_Callback(hObject, eventdata, handles)
% hObject    handle to RunOF (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles = RunOpticalFlowAnalysisButton(handles);
guidata(gcbo,handles);
n=0;


% --- Executes on button press in RunSS.
function RunSS_Callback(hObject, eventdata, handles)
% hObject    handle to RunSS (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
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
n=0;

% --- Executes on button press in ROISpeedAngle.
function ROISpeedAngle_Callback(hObject, eventdata, handles)
% hObject    handle to ROISpeedAngle (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
n=0;



% --- Executes on button press in ViewInstSpeedAngle.
function ViewInstSpeedAngle_Callback(hObject, eventdata, handles)
% hObject    handle to ViewInstSpeedAngle (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles = InstSpDir(handles, eventdata);
if isfield(handles,'ROI')
    if isfield(handles.ROI,'selected')
        if handles.ROI.selected
            plotflag = 1;
            handles = plotInstSpDir(handles, eventdata, plotflag);
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
n = 0;


function handles = InstSpDir(handles,eventdata)
%ViewInstSpeedAngle

methods = {'CLG','HS'};

for method = 1:2
    method = methods{method};
    if isfield(handles,['uv',method,'calculated'])
        if eval(['handles.uv',method,'calculated'])
            handles = get_InstSpDir(handles,eventdata,method);
        end
    end
end

function handles = get_InstSpDir(handles,eventdata,method)
if isfield(handles,['uv',method,'calculated'])
    if eval(['handles.uv',method,'calculated'])
        [~, flag] = getimage(handles.MainAxes); % check if there is an image in the main axes
        if flag
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
                handles = FstartROI_Callback(handles.FstartROI, eventdata, handles);
                handles = FendROI_Callback(handles.FendROI, eventdata, handles);
                FstartROI = max([handles.ROI.Fstart, handles.FstartOFcalculated]);
                FendROI = min([handles.ROI.Fend,handles.FendOFcalculated]);
                set(handles.FstartROI,'string',FstartROI);
                set(handles.FendROI,'string',FendROI);
                
                xy = handles.ROI.xy;
                xyind = sub2ind([handles.dim1, handles.dim2],xy(2,:),xy(1,:));
                nPoints = size(xy,2);
                eval(['handles.InstSpDir',method,' = zeros(nPoints,FendROI-FstartROI);']);
                for idx = FstartROI:FendROI-1
                    eval(['thisframe = handles.uv',method,'(:,:,idx);']);
                    eval(['handles.InstSpDir',method,'(:,idx) = thisframe(xyind);']);
                end
                eval(['handles.FstartROIuv',method,'calculated = FstartROI;']);
                eval(['handles.FendROIuv',method,'calculated = FendROI;']);
                eval(['handles.ROIuv',method,'calculated = xy;']);
            end
        end
    end
end

function handles = plotInstSpDir(handles, eventdata, plotflag)

handles = NbinSpDir_Callback(handles.NbinSpDir, eventdata, handles);
Nbins = handles.HistRose.NbinSpDir;

handles.HistRose.minSp = inf;
handles.HistRose.maxSp = 0;
if isfield(handles,'InstSpDirCLG')
    handles.InstSpDirCLG(isnan(handles.InstSpDirCLG))=0;
    handles.InstSpDirCLG(isinf(handles.InstSpDirCLG))=0;
    
    handles.HistRose.SpCLG = abs(handles.InstSpDirCLG);
    handles.HistRose.minSpCLG = quantile(handles.HistRose.SpCLG(:),0.001);
    handles.HistRose.maxSpCLG = quantile(handles.HistRose.SpCLG(:),0.999);
    handles.HistRose.AngCLG = angle(conj(handles.InstSpDirCLG));
    
    handles.HistRose.minSp = min(handles.HistRose.minSp,handles.HistRose.minSpCLG);
    handles.HistRose.maxSp = max(handles.HistRose.maxSp,handles.HistRose.maxSpCLG);
end
if isfield(handles,'InstSpDirHS')
    handles.InstSpDirHS(isnan(handles.InstSpDirHS))=0;
    handles.InstSpDirHS(isinf(handles.InstSpDirHS))=0;
    
    handles.HistRose.SpHS = abs(handles.InstSpDirHS);
    handles.HistRose.minSpHS = quantile(handles.HistRose.SpHS(:),0.001);
    handles.HistRose.maxSpHS = quantile(handles.HistRose.SpHS(:),0.999);
    handles.HistRose.AngHS = angle(conj(handles.InstSpDirHS));
    
    handles.HistRose.minSp = min(handles.HistRose.minSp,handles.HistRose.minSpHS);
    handles.HistRose.maxSp = max(handles.HistRose.maxSp,handles.HistRose.maxSpHS);
end

% calculate bin centers
binSize = (handles.HistRose.maxSp - handles.HistRose.minSp)/Nbins;
bins = (binSize/2+handles.HistRose.minSp):binSize:handles.HistRose.maxSp;

binsAll = [];
if isfield(handles,'InstSpDirCLG')
    [binsCLG, ~] = hist(handles.HistRose.SpCLG(:),bins);
    binsCLG = 100*binsCLG/length(handles.HistRose.SpCLG(:));
    binsAll = [binsAll;binsCLG];
    handles.InstSpCLG_ele = binsCLG;
    handles.InstSpCLG_bins = bins;
    
    [tCLG, rCLG] = rose(handles.HistRose.AngCLG(:),Nbins);
    rCLG = rCLG/max(rCLG);
    handles.InstDirCLG_rho = rCLG;
    handles.InstDirCLG_theta = tCLG;
end
if isfield(handles,'InstSpDirHS')
    [binsHS, ~] = hist(handles.HistRose.SpHS(:),bins);
    binsHS = 100*binsHS/length(handles.HistRose.SpHS(:));
    binsAll = [binsAll;binsHS];
    handles.InstSpHS_ele = binsHS;
    handles.InstSpHS_bins = bins;
    
    [tHS, rHS] = rose(handles.HistRose.AngHS(:),Nbins);
    rHS = rHS/max(rHS);
    handles.InstDirHS_rho = rHS;
    handles.InstDirHS_theta = tHS;
end

% plotting
if plotflag && (isfield(handles,'InstSpDirCLG') || isfield(handles,'InstSpDirHS'))
    if isfield(handles,'InstSpDirFig')
        if ishandle(handles.InstSpDirFig)
            figure(handles.InstSpDirFig);
        else
            handles.InstSpDirFig = figure(102);
            set(handles.InstSpDirFig,'visible','off');
            set(handles.InstSpDirFig,'numbertitle','off','name','Instantaneous Speeds and Angles histogram for selected ROI','units','inches','position',[1 3.5 6 2.5]);
            set(handles.InstSpDirFig,'visible','on');
        end
    else
        handles.InstSpDirFig = figure(102);
        set(handles.InstSpDirFig,'visible','off');
        set(handles.InstSpDirFig,'numbertitle','off','name','Instantaneous Speeds and Angles histogram for selected ROI','units','inches','position',[1 3.5 6 2.5]);
        set(handles.InstSpDirFig,'visible','on');
    end
    c1 = [0.08,0.17,0.55];
    c2 = [0.55,0,0];
    figure(handles.InstSpDirFig);
    subplot(121); cla;
    bar(bins,binsAll');
    maxY = max(binsAll(:));
    xlim([min(bins)-binSize, max(bins)+binSize]);
    ylim([0 maxY*1.1]);
    set(gca,'box','off','TickDir','out');
    xlabel(sprintf('S_i (p/f)')); ylabel('Percentage');
    
    figure(handles.InstSpDirFig);
    subplot(122); cla;
    if isfield(handles,'InstSpDirCLG')
        h = polar(tCLG,rCLG);
        set(h,'color',c1);
        hold on;
        if isfield(handles,'InstSpDirHS')
            h = polar(tHS,rHS); 
            set(h,'color',c2);
            legend('CLG','HS','Location','best')
            hold off;
        else
            legend('CLG','Location','best')
        end
    else
        if isfield(handles,'InstSpDirHS')
            h = polar(tHS,rHS); 
            set(h,'color',c1);
            legend('HS','Location','best')
            hold off;
        end
    end
end


% --- Executes on button press in pushbutton11.
function pushbutton11_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton11 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in checkbox5.
function checkbox5_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox5


% --- Executes on button press in LoadImgSeq.
function LoadImgSeq_Callback(hObject, eventdata, handles)
% hObject    handle to LoadImgSeq (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
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
% hObject    handle to LoadMask (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
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
% hObject    handle to MaxLagTS (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of MaxLagTS as text
%        str2double(get(hObject,'String')) returns contents of MaxLagTS as a double


% --- Executes during object creation, after setting all properties.
function MaxLagTS_CreateFcn(hObject, eventdata, handles)
% hObject    handle to MaxLagTS (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in MaskChk.
function MaskChk_Callback(hObject, eventdata, handles)
% hObject    handle to MaskChk (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of MaskChk
currentFrame = round(get(handles.FramesSlider,'Value'));
handles = loadFrame(handles,eventdata,currentFrame);
guidata(gcbo,handles);
n=0;

% --- Executes on button press in IntensityTime.
function handles = IntensityTime_Callback(hObject, eventdata, handles)
% hObject    handle to IntensityTime (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
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
    
    handles = FstartROI_Callback(handles.FstartROI, eventdata, handles);
    handles = FendROI_Callback(handles.FendROI, eventdata, handles);
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
% hObject    handle to ViewTempSpeedTrajLength (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles = TempSpTrajLen(handles,eventdata);
if isfield(handles,'ROI')
    if isfield(handles.ROI,'selected')
        if handles.ROI.selected
            plotflag = 1;
            handles = plotTempSpTrajLen(handles, eventdata, plotflag);
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

function handles = TempSpTrajLen(handles,eventdata)
methods = {'CLG','HS'};
for method = 1:2
    method = methods{method};
    if isfield(handles,['uv',method,'calculated'])
        if eval(['handles.uv',method,'calculated'])
            handles = get_TempSpTrajLen(handles,eventdata,method);
        end
    end
end

function handles = get_TempSpTrajLen(handles,eventdata,method)
if isfield(handles,['uv',method,'calculated'])
    if eval(['handles.uv',method,'calculated'])
        [~, flag] = getimage(handles.MainAxes); % check if there is an image in the main axes
        if flag
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
                handles = FstartROI_Callback(handles.FstartROI, eventdata, handles);
                handles = FendROI_Callback(handles.FendROI, eventdata, handles);
                FstartROI = max([handles.ROI.Fstart, handles.FstartOFcalculated]);
                FendROI = min([handles.ROI.Fend,handles.FendOFcalculated]);
                set(handles.FstartROI,'string',FstartROI);
                set(handles.FendROI,'string',FendROI);
                
                xy = handles.ROI.xy;
                sx = xy(1,:);
                sy = xy(2,:);
                startFrame = FstartROI - handles.FstartOFcalculated +1;
                endFrame = FendROI - handles.FstartOFcalculated +1;
                offsetFrame = handles.FstartOFcalculated;
                if startFrame < endFrame
                    output = Velocity_Profile(eval(['handles.uv',method]),startFrame,endFrame,offsetFrame,sx,sy);
                    
                    eval(['handles.TempSpTrajLen',method,' = output;']);
                    eval(['handles.FstartROITemp',method,'calculated = FstartROI;']);
                    eval(['handles.FendROITemp',method,'calculated = FendROI;']);
                    eval(['handles.ROITemp',method,'calculated = xy;']);
                end
            end
        end
    end
end

function handles = plotTempSpTrajLen(handles, eventdata, plotflag)

handles = NbinSpTraj_Callback(handles.NbinSpTraj, eventdata, handles);
Nbins = handles.TrajHist.NbinSpTraj;

handles.TrajHist.minSp = inf;
handles.TrajHist.maxSp = 0;
handles.TrajHist.minLen = inf;
handles.TrajHist.maxLen = 0;
handles.TrajHist.minMaxSp = inf;
handles.TrajHist.maxMaxSp = 0;
if isfield(handles,'TempSpTrajLenCLG')
    AllSpCLG = zeros(length(handles.TempSpTrajLenCLG.speed{1}),length(handles.TempSpTrajLenCLG.speed));
    totalDistCLG = zeros(1,length(handles.TempSpTrajLenCLG.speed));
    for point = 1:length(handles.TempSpTrajLenCLG.speed)
        AllSpCLG(:,point) = handles.TempSpTrajLenCLG.speed{point};
        
        str = handles.TempSpTrajLenCLG.str{point};
        xs = str(:,1); ys = str(:,2);
        vects = xs + 1i*ys;
        diffVects = diff(vects);
        distances = abs(diffVects);
        totalDistCLG(point) = sum(distances);
    end
    
    tSpCLG = handles.TempSpTrajLenCLG.tspeed{point};
    AllSpCLG(isnan(AllSpCLG))=0;
    AllSpCLG(isinf(AllSpCLG))=0;
    totalDistCLG(isnan(totalDistCLG))=0;
    totalDistCLG(isinf(totalDistCLG))=0;
    AllMaxSpCLG = max(AllSpCLG);
    
    handles.TrajHist.minSpCLG = quantile(AllSpCLG(:),0.001);
    handles.TrajHist.maxSpCLG = quantile(AllSpCLG(:),0.999);
    handles.TrajHist.minSp = min(handles.TrajHist.minSp,handles.TrajHist.minSpCLG);
    handles.TrajHist.maxSp = max(handles.TrajHist.maxSp,handles.TrajHist.maxSpCLG);
    
    handles.TrajHist.minLenCLG = quantile(totalDistCLG(:),0.001);
    handles.TrajHist.maxLenCLG = quantile(totalDistCLG(:),0.999);
    handles.TrajHist.minLen = min(handles.TrajHist.minLen,handles.TrajHist.minLenCLG);
    handles.TrajHist.maxLen = max(handles.TrajHist.maxLen,handles.TrajHist.maxLenCLG);
    
    handles.TrajHist.minMaxSpCLG = quantile(AllMaxSpCLG,0.001);
    handles.TrajHist.maxMaxSpCLG = quantile(AllMaxSpCLG,0.999);
    handles.TrajHist.minMaxSp = min(handles.TrajHist.minMaxSp,handles.TrajHist.minMaxSpCLG);
    handles.TrajHist.maxMaxSp = max(handles.TrajHist.maxMaxSp,handles.TrajHist.maxMaxSpCLG);
end

if isfield(handles,'TempSpTrajLenHS')
    AllSpHS = zeros(length(handles.TempSpTrajLenHS.speed{1}),length(handles.TempSpTrajLenHS.speed));
    totalDistHS = zeros(1,length(handles.TempSpTrajLenHS.speed));
    for point = 1:length(handles.TempSpTrajLenHS.speed)
        AllSpHS(:,point) = handles.TempSpTrajLenHS.speed{point};
        
        str = handles.TempSpTrajLenHS.str{point};
        xs = str(:,1); ys = str(:,2);
        vects = xs + 1i*ys;
        diffVects = diff(vects);
        distances = abs(diffVects);
        totalDistHS(point) = sum(distances);
    end
    AllMaxSpHS = max(AllSpHS);
    tSpHS = handles.TempSpTrajLenHS.tspeed{point};
    AllSpHS(isnan(AllSpHS))=0;
    AllSpHS(isinf(AllSpHS))=0;
    totalDistHS(isnan(totalDistHS))=0;
    totalDistHS(isinf(totalDistHS))=0;
    
    handles.TrajHist.minSpHS = quantile(AllSpHS(:),0.001);
    handles.TrajHist.maxSpHS = quantile(AllSpHS(:),0.999);
    handles.TrajHist.minSp = min(handles.TrajHist.minSp,handles.TrajHist.minSpHS);
    handles.TrajHist.maxSp = max(handles.TrajHist.maxSp,handles.TrajHist.maxSpHS);
    
    handles.TrajHist.minLenHS = quantile(totalDistHS(:),0.001);
    handles.TrajHist.maxLenHS = quantile(totalDistHS(:),0.999);
    handles.TrajHist.minLen = min(handles.TrajHist.minLen,handles.TrajHist.minLenHS);
    handles.TrajHist.maxLen = max(handles.TrajHist.maxLen,handles.TrajHist.maxLenHS);
    
    handles.TrajHist.minMaxSpHS = quantile(AllMaxSpHS,0.001);
    handles.TrajHist.maxMaxSpHS = quantile(AllMaxSpHS,0.999);
    handles.TrajHist.minMaxSp = min(handles.TrajHist.minMaxSp,handles.TrajHist.minMaxSpHS);
    handles.TrajHist.maxMaxSp = max(handles.TrajHist.maxMaxSp,handles.TrajHist.maxMaxSpHS);
end

% calculate bin centers
binSizeSp = (handles.TrajHist.maxSp - handles.TrajHist.minSp)/Nbins;
binsSp = (binSizeSp/2+handles.TrajHist.minSp):binSizeSp:handles.TrajHist.maxSp;

binSizeLen = (handles.TrajHist.maxLen - handles.TrajHist.minLen)/Nbins;
binsLen = (binSizeLen/2+handles.TrajHist.minLen):binSizeLen:handles.TrajHist.maxLen;

binSizeMaxSp = (handles.TrajHist.maxMaxSp - handles.TrajHist.minMaxSp)/Nbins;
binsMaxSp = (binSizeMaxSp/2+handles.TrajHist.minMaxSp):binSizeMaxSp:handles.TrajHist.maxMaxSp;

binsAllSp = [];
binsAllLen = [];
binsAllMaxSp = [];
if isfield(handles,'TempSpTrajLenCLG')
    [binsSpCLG, ~] = hist(AllSpCLG(:),binsSp);
    binsSpCLG = 100*binsSpCLG/length(AllSpCLG(:));
    binsAllSp = [binsAllSp;binsSpCLG];
    handles.TempSpCLG_ele = binsSpCLG;
    handles.TempSpCLG_bins = binsSp;
    
    [binsLenCLG, ~] = hist(totalDistCLG,binsLen);
    binsLenCLG = 100*binsLenCLG/length(totalDistCLG);
    binsAllLen = [binsAllLen;binsLenCLG];
    handles.TempLenCLG_ele = binsLenCLG;
    handles.TempLenCLG_bins = binsLen;
    
    [binsMaxCLG, ~] = hist(AllMaxSpCLG(:),binsMaxSp);
    binsMaxCLG = 100*binsMaxCLG/length(AllMaxSpCLG(:));
    binsAllMaxSp = [binsAllMaxSp;binsMaxCLG];
    handles.MaxTempSpCLG_ele = binsMaxCLG;
    handles.MaxTempSpCLG_bins = binsMaxSp;
end

if isfield(handles,'TempSpTrajLenHS')
    [binsSpHS, ~] = hist(AllSpHS(:),binsSp);
    binsSpHS = 100*binsSpHS/length(AllSpHS(:));
    binsAllSp = [binsAllSp;binsSpHS];
    handles.TempSpHS_ele = binsSpHS;
    handles.TempSpHS_bins = binsSp;
    
    [binsLenHS, ~] = hist(totalDistHS,binsLen);
    binsLenHS = 100*binsLenHS/length(totalDistHS);
    binsAllLen = [binsAllLen;binsLenHS];
    handles.TempLenHS_ele = binsLenHS;
    handles.TempLenHS_bins = binsLen;
    
    [binsMaxHS, ~] = hist(AllMaxSpHS(:),binsMaxSp);
    binsMaxHS = 100*binsMaxHS/length(AllMaxSpHS(:));
    binsAllMaxSp = [binsAllMaxSp;binsMaxHS];
    handles.MaxTempSpHS_ele = binsMaxHS;
    handles.MaxTempSpHS_bins = binsMaxSp;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plotting
if plotflag && (isfield(handles,'TempSpTrajLenCLG') || isfield(handles,'TempSpTrajLenHS'))
    if isfield(handles,'TempSptrajLenFig')
        if ishandle(handles.TempSptrajLenFig)
            figure(handles.TempSptrajLenFig);
        else
            handles.TempSptrajLenFig = figure(103);
            set(handles.TempSptrajLenFig,'visible','off');
            set(handles.TempSptrajLenFig,'numbertitle','off','name','Temporal Speeds, profile and histogram, and Trajectory Length histogram for selected ROI','units','inches','position',[1 1 12 2.5]);
            set(handles.TempSptrajLenFig,'visible','on');
        end
    else
        handles.TempSptrajLenFig = figure(103);
        set(handles.TempSptrajLenFig,'visible','off');
        set(handles.TempSptrajLenFig,'numbertitle','off','name','Temporal Speeds, profile and histogram, and Trajectory Length histogram for selected ROI','units','inches','position',[1 1 12 2.5]);
        set(handles.TempSptrajLenFig,'visible','on');
    end
    clf
    % temporal profile
    c1 = [0.08,0.17,0.55];
    c2 = [0.55,0,0];
    figure(handles.TempSptrajLenFig);
    subplot(141); cla;
    if isfield(handles,'TempSpTrajLenCLG')
        LegendCLG = plot(tSpCLG,AllSpCLG,'color',c1);hold on;
        if isfield(handles,'TempSpTrajLenHS')
            LegendHS = plot(tSpHS,AllSpHS,'color',c2); 
            legend([LegendCLG(1),LegendHS(1)],'CLG','HS','Location','best')
            hold off;
        else
            legend('CLG','Location','best')
        end
    else
        if isfield(handles,'TempSpTrajLenHS')
            plot(tSpHS,AllSpHS,'color',c1); 
            legend('HS','Location','best')
            hold off;
        end
    end
    ylabel('S_t_e (p/f)');
    xlabel('Frame Number');
    
    % teporal speed histogram
    figure(handles.TempSptrajLenFig);
    subplot(142); cla;
    
    bar(binsSp,binsAllSp');
    maxY = max(binsAllSp(:));
    xlim([min(binsSp)-binSizeSp, max(binsSp)+binSizeSp]);
    ylim([0 maxY*1.1]);
    set(gca,'box','off','TickDir','out');
    xlabel(sprintf('S_t_e (p/f)')); ylabel('Percentage');
    
    % Maximum teporal speed histogram
    figure(handles.TempSptrajLenFig);
    subplot(143); cla;
    
    bar(binsMaxSp,binsAllMaxSp');
    maxY = max(binsAllMaxSp(:));
    xlim([min(binsMaxSp)-binSizeMaxSp, max(binsMaxSp)+binSizeMaxSp]);
    ylim([0 maxY*1.1]);
    set(gca,'box','off','TickDir','out');
    xlabel(sprintf('Max S_t_e (p/f)')); ylabel('Percentage');
    
    % Trajectory Distance histogram
    figure(handles.TempSptrajLenFig);
    subplot(144); cla;
    
    bar(binsLen,binsAllLen');
    maxY = max(binsAllLen(:));
    xlim([min(binsLen)-binSizeLen, max(binsLen)+binSizeLen]);
    ylim([0 maxY*1.1]);
    set(gca,'box','off','TickDir','out');
    xlabel(sprintf('Trajectory Length (pixels)')); ylabel('Percentage');
end


% --- Executes on button press in SaveSpeedAngle.
function SaveSpeedAngle_Callback(hObject, eventdata, handles)
% hObject    handle to SaveSpeedAngle (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Save ROI_InstSpAng_TempSpTraj
handles = FstartROI_Callback(handles.FstartROI, eventdata, handles);
handles = FendROI_Callback(handles.FendROI, eventdata, handles);
if isfield(handles,'FstartOFcalculated')
    FstartROI = max([handles.ROI.Fstart, handles.FstartOFcalculated]);
else
    FstartROI = handles.ROI.Fstart;
end
if isfield(handles,'FendOFcalculated')
    FendROI = min([handles.ROI.Fend,handles.FendOFcalculated]);
else
    FendROI = handles.ROI.Fend;
end
set(handles.FstartROI,'string',FstartROI);
set(handles.FendROI,'string',FendROI);

% CLG InstSpAng
if isfield(handles,'InstSpAngcalculated')
    if isfield(handles,'ROI')
        if isfield(handles.ROI,'selected')
            if handles.ROI.selected && isfield(handles,'FstartROIuvCLGcalculated') &&...
                    isfield(handles,'FendROIuvCLGcalculated') && isfield(handles,'ROIuvCLGcalculated')
                if (handles.ROI.Fstart == handles.FstartROIuvCLGcalculated) &&...
                        (handles.ROI.Fend == handles.FendROIuvCLGcalculated) &&...
                        (size(handles.ROI.xy,2) == size(handles.ROIuvCLGcalculated,2))
                    if handles.ROI.xy == handles.ROIuvCLGcalculated
                        saveInstSpAngCLG = 1;
                    else
                        saveInstSpAngCLG = 0; % go and calc inst for CLG
                    end
                else
                    saveInstSpAngCLG = 0; % go and calc inst for CLG
                end
                saveInstSpAngCLG = 1;
            else
                saveInstSpAngCLG = 0; % go and calc inst for CLG
            end
        else
            saveInstSpAngCLG = 0; % go and calc inst for CLG
        end
    end
else
    if isfield(handles,'uvCLGcalculated')
        if  handles.uvCLGcalculated
            %             if ~isfield(handles,'ROI')
            %                 handles = SelectROI(handles.SelectROI, eventdata, handles);
            %             end
            saveInstSpAngCLG = 0; % go and calc inst for CLG
        end
    end
end

% HS InstSpAng
if isfield(handles,'InstSpAngcalculated')
    if isfield(handles,'ROI')
        if isfield(handles.ROI,'selected')
            if handles.ROI.selected && isfield(handles,'FstartROIuvHScalculated') &&...
                    isfield(handles,'FendROIuvHScalculated') && isfield(handles,'ROIuvHScalculated')
                if (handles.ROI.Fstart == handles.FstartROIuvHScalculated) &&...
                        (handles.ROI.Fend == handles.FendROIuvHScalculated) &&...
                        (size(handles.ROI.xy,2) == size(handles.ROIuvCLGcalculated,2))
                    if handles.ROI.xy == handles.ROIuvCLGcalculated
                        saveInstSpAngHS = 1;
                    else
                        saveInstSpAngHS = 0; % go and calc inst for HS
                    end
                else
                    saveInstSpAngHS = 0; % go and calc inst for HS
                end
                saveInstSpAngHS = 1;
            else
                saveInstSpAngHS = 0; % go and calc inst for HS
            end
        else
            saveInstSpAngHS = 0; % go and calc inst for HS
        end
    end
else
    if isfield(handles,'uvHScalculated')
        if  handles.uvHScalculated
            %             if ~isfield(handles,'ROI')
            %                 handles = SelectROI(handles.SelectROI, eventdata, handles);
            %             end
            saveInstSpAngHS = 0; % go and calc inst for HS
        end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CLG TempSpTraj
if isfield(handles,'TempSpTrajLengthcalculated')
    if isfield(handles,'ROI')
        if isfield(handles.ROI,'selected')
            if handles.ROI.selected && isfield(handles,'FstartROITempCLGcalculated') &&...
                    isfield(handles,'FendROITempCLGcalculated') && isfield(handles,'ROITempCLGcalculated')
                if (handles.ROI.Fstart == handles.FstartROITempCLGcalculated) &&...
                        (handles.ROI.Fend == handles.FendROITempCLGcalculated) &&...
                        (size(handles.ROI.xy,2) == size(handles.ROITempCLGcalculated,2))
                    if handles.ROI.xy == handles.ROITempCLGcalculated
                        saveTempSpTrajCLG = 1;
                    else
                        saveTempSpTrajCLG = 0; % go and calc temp for CLG
                    end
                else
                    saveTempSpTrajCLG = 0; % go and calc temp for CLG
                end
                saveTempSpTrajCLG = 1;
            else
                saveTempSpTrajCLG = 0; % go and calc temp for CLG
            end
        else
            saveTempSpTrajCLG = 0; % go and calc temp for CLG
        end
    end
else
    if isfield(handles,'uvCLGcalculated')
        if  handles.uvCLGcalculated
            %             if ~isfield(handles,'ROI')
            %                 handles = SelectROI(handles.SelectROI, eventdata, handles);
            %             end
            saveTempSpTrajCLG = 0; % go and calc temp for CLG
        end
    end
end

% HS TempSpTraj
if isfield(handles,'TempSpTrajLengthcalculated')
    if isfield(handles,'ROI')
        if isfield(handles.ROI,'selected')
            if handles.ROI.selected  && isfield(handles,'FstartROITempCLGcalculated') &&...
                    isfield(handles,'FendROITempCLGcalculated') && isfield(handles,'ROITempCLGcalculated')
                if (handles.ROI.Fstart == handles.FstartROITempHScalculated) &&...
                        (handles.ROI.Fend == handles.FendROITempHScalculated) &&...
                        (size(handles.ROI.xy,2) == size(handles.ROITempHScalculated,2))
                    if handles.ROI.xy == handles.ROITempHScalculated
                        saveTempSpTrajHS = 1;
                    else
                        saveTempSpTrajHS = 0; % go and calc temp for HS
                    end
                else
                    saveTempSpTrajHS = 0; % go and calc temp for HS
                end
                saveTempSpTrajHS = 1;
            else
                saveTempSpTrajHS = 0; % go and calc temp for HS
            end
        else
            saveTempSpTrajHS = 0; % go and calc temp for HS
        end
    end
else
    if isfield(handles,'uvHScalculated')
        if  handles.uvHScalculated
            %             if ~isfield(handles,'ROI')
            %                 handles = SelectROI(handles.SelectROI, eventdata, handles);
            %             end
            saveTempSpTrajHS = 0; % go and calc temp for HS
        end
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% go and calc inst for CLG & HS
if isfield(handles,'uvCLGcalculated') || isfield(handles,'uvHScalculated')
    if  handles.uvCLGcalculated || handles.uvHScalculated
        if saveInstSpAngCLG == 0 || saveInstSpAngHS == 0
            handles = InstSpDir(handles, eventdata);
            if handles.ROI.selected
                plotflag = 0;
                handles = plotInstSpDir(handles, eventdata, plotflag);
                handles.InstSpAngcalculated = 1;
            end
        end
    end
end

% go and calc temp for CLG & HS
if isfield(handles,'uvCLGcalculated') || isfield(handles,'uvHScalculated')
    if  handles.uvCLGcalculated || handles.uvHScalculated
        if saveTempSpTrajCLG == 0 || saveTempSpTrajHS == 0
            handles = TempSpTrajLen(handles,eventdata);
            if handles.ROI.selected
                plotflag = 1;
                handles = plotTempSpTrajLen(handles, eventdata, plotflag);
                handles.TempSpTrajLengthcalculated = 1;
            end
        end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ~isfield(handles.ROI,'selected')
    handles.ROI.selected = 0;
end
if handles.ROI.selected
    handles = ROISaveNameSpeedAngle_Callback(handles.ROISaveNameSpeedAngle, eventdata, handles);
    FileName_save_SpAngTraj = handles.FileName_save_SpAngTraj;
    
    if isfield(handles,'PathName')
        handles.SavePathName = handles.PathName;
    else
        handles.SavePathName = pwd; % current folder
    end
    % save for CLG
    if isfield(handles,'uvCLGcalculated')
        if  handles.uvCLGcalculated
            if ~exist(handles.SavePathName,'dir')
                mkdir(handles.SavePathName);
            end
            SaveFullFileName = [handles.SavePathName,'\',FileName_save_SpAngTraj,'.mat'];
            if exist(SaveFullFileName)
                load(SaveFullFileName);
            end
            mFileROI_InstSpAng_TempSpTraj = matfile(SaveFullFileName,'Writable',true);
            
            InstSpAng.CLG.uv = handles.InstSpDirCLG;
            InstSpAng.CLG.SpHist.ele = handles.InstSpCLG_ele;
            InstSpAng.CLG.SpHist.bins = handles.InstSpCLG_bins;
            InstSpAng.CLG.AngHist.rho = handles.InstDirCLG_rho;
            InstSpAng.CLG.AngHist.theta = handles.InstDirCLG_theta;
            
            mFileROI_InstSpAng_TempSpTraj.InstSpAng = InstSpAng;
            
            TempSpTraj.CLG.StrLines = handles.TempSpTrajLenCLG;
            TempSpTraj.CLG.TempSpHist.ele = handles.TempSpCLG_ele;
            TempSpTraj.CLG.TempSpHist.bins = handles.TempSpCLG_bins;
            TempSpTraj.CLG.TempLenHist.ele = handles.TempLenCLG_ele;
            TempSpTraj.CLG.TempLenHist.bins = handles.TempLenCLG_bins;
            TempSpTraj.CLG.MaxTempSpHist.ele = handles.MaxTempSpCLG_ele;
            TempSpTraj.CLG.MaxTempSpHist.bins = handles.MaxTempSpCLG_bins;
            
            mFileROI_InstSpAng_TempSpTraj.TempSpTraj = TempSpTraj;
        end
    end
    
    % save for HS
    if isfield(handles,'uvHScalculated')
        if  handles.uvHScalculated
            if ~exist(handles.SavePathName,'dir')
                mkdir(handles.SavePathName);
            end
            SaveFullFileName = [handles.SavePathName,'\',FileName_save_SpAngTraj,'.mat'];
            if exist(SaveFullFileName)
                load(SaveFullFileName);
            end
            mFileROI_InstSpAng_TempSpTraj = matfile(SaveFullFileName,'Writable',true);
            
            InstSpAng.HS.uv = handles.InstSpDirHS;
            InstSpAng.HS.SpHist.ele = handles.InstSpHS_ele;
            InstSpAng.HS.SpHist.bins = handles.InstSpHS_bins;
            InstSpAng.HS.AngHist.rho = handles.InstDirHS_rho;
            InstSpAng.HS.AngHist.theta = handles.InstDirHS_theta;
            
            mFileROI_InstSpAng_TempSpTraj.InstSpAng = InstSpAng;
            
            TempSpTraj.HS.StrLines = handles.TempSpTrajLenHS;
            TempSpTraj.HS.TempSpHist.ele = handles.TempSpHS_ele;
            TempSpTraj.HS.TempSpHist.bins = handles.TempSpHS_bins;
            TempSpTraj.HS.TempLenHist.ele = handles.TempLenHS_ele;
            TempSpTraj.HS.TempLenHist.bins = handles.TempLenHS_bins;
            TempSpTraj.HS.MaxTempSpHist.ele = handles.MaxTempSpHS_ele;
            TempSpTraj.HS.MaxTempSpHist.bins = handles.MaxTempSpHS_bins;
            
            mFileROI_InstSpAng_TempSpTraj.TempSpTraj = TempSpTraj;
        end
    end
    
end

guidata(gcbo,handles);
n=0;


% --- Executes on button press in ViewTraj.
function ViewTraj_Callback(hObject, eventdata, handles)
% hObject    handle to ViewTraj (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)



function handles = ROISaveNameSpeedAngle_Callback(hObject, eventdata, handles)
% hObject    handle to ROISaveNameSpeedAngle (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of ROISaveNameSpeedAngle as text
%        str2double(get(hObject,'String')) returns contents of ROISaveNameSpeedAngle as a double
handles.FileName_save_SpAngTraj = get(hObject,'string');
handles.FileName_save_SpAngTraj = strtrim(handles.FileName_save_SpAngTraj);
if strcmp(handles.FileName_save_SpAngTraj,'')
    handles.FileName_save_SpAngTraj = 'ROI_InstSpAng_TempSpTraj';
    set(hObject,'string',handles.FileName_save_SpAngTraj);
end
guidata(gcbo,handles);



% --- Executes during object creation, after setting all properties.
function ROISaveNameSpeedAngle_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ROISaveNameSpeedAngle (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in SSChk.
function SSChk_Callback(hObject, eventdata, handles)
% hObject    handle to SSChk (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of SSChk
currentFrame = round(get(handles.FramesSlider,'Value'));
handles = loadFrame(handles,eventdata,currentFrame);
guidata(gcbo,handles);

% --- Executes on button press in SaveSS.
function SaveSS_Callback(hObject, eventdata, handles)
% hObject    handle to SaveSS (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of SaveSS


% --- Executes on button press in SaveOF.
function SaveOF_Callback(hObject, eventdata, handles)
% hObject    handle to SaveOF (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of SaveOF


% --- Executes on button press in HSsave.
function HSsave_Callback(hObject, eventdata, handles)
% hObject    handle to HSsave (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of HSsave


% --- Executes on button press in CLGsave.
function CLGsave_Callback(hObject, eventdata, handles)
% hObject    handle to CLGsave (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of CLGsave


% --- Executes on button press in TSsave.
function TSsave_Callback(hObject, eventdata, handles)
% hObject    handle to TSsave (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of TSsave


% --- Executes on button press in OFrawFile.
function OFrawFile_Callback(hObject, eventdata, handles)
% hObject    handle to OFrawFile (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of OFrawFile


% --- Executes on button press in OFmatFile.
function OFmatFile_Callback(hObject, eventdata, handles)
% hObject    handle to OFmatFile (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of OFmatFile



function FstartOF_Callback(hObject, eventdata, handles)
% hObject    handle to FstartOF (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of FstartOF as text
%        str2double(get(hObject,'String')) returns contents of FstartOF as a double


% --- Executes during object creation, after setting all properties.
function FstartOF_CreateFcn(hObject, eventdata, handles)
% hObject    handle to FstartOF (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function FendOF_Callback(hObject, eventdata, handles)
% hObject    handle to FendOF (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of FendOF as text
%        str2double(get(hObject,'String')) returns contents of FendOF as a double


% --- Executes during object creation, after setting all properties.
function FendOF_CreateFcn(hObject, eventdata, handles)
% hObject    handle to FendOF (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function handles = MinDispVal_Callback(hObject, eventdata, handles)
% hObject    handle to MinDispVal (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of MinDispVal as text
%        str2double(get(hObject,'String')) returns contents of MinDispVal as a double
MinDispVal_default = 0;
handles.MinVal = Str2NumFromHandle(hObject,MinDispVal_default);
currentFrame = round(get(handles.FramesSlider,'Value'));
handles = loadFrame(handles,eventdata,currentFrame);
guidata(gcbo,handles);


% --- Executes during object creation, after setting all properties.
function MinDispVal_CreateFcn(hObject, eventdata, handles)
% hObject    handle to MinDispVal (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function handles = MaxDispVal_Callback(hObject, eventdata, handles)
% hObject    handle to MaxDispVal (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of MaxDispVal as text
%        str2double(get(hObject,'String')) returns contents of MaxDispVal as a double
MaxDispVal_default = 0;
handles.MaxVal = Str2NumFromHandle(hObject,MaxDispVal_default);
currentFrame = round(get(handles.FramesSlider,'Value'));
handles = loadFrame(handles,eventdata,currentFrame);
guidata(gcbo,handles);


% --- Executes during object creation, after setting all properties.
function MaxDispVal_CreateFcn(hObject, eventdata, handles)
% hObject    handle to MaxDispVal (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in checkbox17.
function checkbox17_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox17 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox17


% --- Executes on button press in SelectROI.
function handles = SelectROI_Callback(hObject, eventdata, handles)
% hObject    handle to SelectROI (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles = SelectROI(hObject, eventdata, handles);
guidata(gcbo,handles);


% --- Executes on button press in SelectPoints.
function handles = SelectPoints_Callback(hObject, eventdata, handles)
% hObject    handle to SelectPoints (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
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

function handles = FstartROI_Callback(hObject, eventdata, handles)
% hObject    handle to FstartROI (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of FstartROI as text
%        str2double(get(hObject,'String')) returns contents of FstartROI as a double
FstartROI_default = 1;
FstartROI = Str2NumFromHandle(handles.FstartROI,FstartROI_default);
FstartROI = round(FstartROI);
set(handles.FstartROI,'string',FstartROI);
handles.ROI.Fstart = FstartROI;
guidata(gcbo,handles);


% --- Executes during object creation, after setting all properties.
function FstartROI_CreateFcn(hObject, eventdata, handles)
% hObject    handle to FstartROI (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function handles = FendROI_Callback(hObject, eventdata, handles)
% hObject    handle to FendROI (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of FendROI as text
%        str2double(get(hObject,'String')) returns contents of FendROI as a double
if isfield(handles,'nFrames')
    FendROI_default = handles.nFrames;
else
    FendROI_default = 1;
end
FendROI = Str2NumFromHandle(handles.FendROI,FendROI_default);
FendROI = round(FendROI);
set(handles.FendROI,'string',FendROI);
handles.ROI.Fend = FendROI;
guidata(gcbo,handles);

% --- Executes during object creation, after setting all properties.
function FendROI_CreateFcn(hObject, eventdata, handles)
% hObject    handle to FendROI (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in checkbox16.
function checkbox16_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox16 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox16


% --- Executes on button press in Trajectory.
function Trajectory_Callback(hObject, eventdata, handles)
% hObject    handle to Trajectory (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of Trajectory
if get(hObject,'value')
    handles.TrajectoryInfo.currVal = 1;
    if get(handles.CLGVisOF,'value'); methodTraj = 'CLG'; end
    if get(handles.HSVisOF,'value'); methodTraj = 'HS'; end
    
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
                handles = FstartROI_Callback(handles.FstartROI, eventdata, handles);
                handles = FendROI_Callback(handles.FendROI, eventdata, handles);
                FstartROI = max([handles.ROI.Fstart, handles.FstartOFcalculated]);
                FendROI = min([handles.ROI.Fend,handles.FendOFcalculated]);
                set(handles.FstartROI,'string',FstartROI);
                set(handles.FendROI,'string',FendROI);
                
                xy = handles.ROI.xy;
                sx = xy(1,:);
                sy = xy(2,:);
                startFrame = FstartROI - handles.FstartOFcalculated +1;
                endFrame = FendROI - handles.FstartOFcalculated +1;
                offsetFrame = handles.FstartOFcalculated;
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


function handles = NbinSpDir_Callback(hObject, eventdata, handles)
% hObject    handle to NbinSpDir (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of NbinSpDir as text
%        str2double(get(hObject,'String')) returns contents of NbinSpDir as a double
NbinSpDir_default = 20;
NbinSpDir = Str2NumFromHandle(handles.NbinSpDir,NbinSpDir_default);
NbinSpDir = round(NbinSpDir);
set(handles.NbinSpDir,'string',NbinSpDir);
handles.HistRose.NbinSpDir = NbinSpDir;
guidata(gcbo,handles);


% --- Executes during object creation, after setting all properties.
function NbinSpDir_CreateFcn(hObject, eventdata, handles)
% hObject    handle to NbinSpDir (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function handles = NbinSpTraj_Callback(hObject, eventdata, handles)
% hObject    handle to NbinSpTraj (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of NbinSpTraj as text
%        str2double(get(hObject,'String')) returns contents of NbinSpTraj as a double
NbinSpTraj_default = 20;
NbinSpTraj = Str2NumFromHandle(handles.NbinSpTraj,NbinSpTraj_default);
NbinSpTraj = round(NbinSpTraj);
set(handles.NbinSpTraj,'string',NbinSpTraj);
handles.TrajHist.NbinSpTraj = NbinSpTraj;
guidata(gcbo,handles);


% --- Executes during object creation, after setting all properties.
function NbinSpTraj_CreateFcn(hObject, eventdata, handles)
% hObject    handle to NbinSpTraj (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in CLGVisSS.
function CLGVisSS_Callback(hObject, eventdata, handles)
% hObject    handle to CLGVisSS (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of CLGVisSS
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
% hObject    handle to HSVisSS (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of HSVisSS
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
% hObject    handle to CLGVisTraj (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of CLGVisTraj
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
% hObject    handle to HSVisTraj (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of HSVisTraj
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
% hObject    handle to CLGVisOF (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of CLGVisOF
if get(hObject,'value')
    set(handles.HSVisOF,'value',0);
else
    set(handles.HSVisOF,'value',1);
end
VectorFieldVisChk_Callback(handles.VectorFieldVisChk, eventdata, handles);
guidata(gcbo,handles);

% --- Executes on button press in HSVisOF.
function HSVisOF_Callback(hObject, eventdata, handles)
% hObject    handle to HSVisOF (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of HSVisOF
if get(hObject,'value')
    set(handles.CLGVisOF,'value',0);
else
    set(handles.CLGVisOF,'value',1);
end
VectorFieldVisChk_Callback(handles.VectorFieldVisChk, eventdata, handles);
guidata(gcbo,handles);


% --- Executes on button press in runCLG.
function runCLG_Callback(hObject, eventdata, handles)
% hObject    handle to runCLG (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of runCLG


% --- Executes on button press in runHS.
function runHS_Callback(hObject, eventdata, handles)
% hObject    handle to runHS (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of runHS


% --- Executes on button press in runTS.
function runTS_Callback(hObject, eventdata, handles)
% hObject    handle to runTS (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of runTS


% --- Executes on button press in SpiralChk.
function SpiralChk_Callback(hObject, eventdata, handles)
% hObject    handle to SpiralChk (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of SpiralChk
currentFrame = round(get(handles.FramesSlider,'Value'));
handles = loadFrame(handles,eventdata,currentFrame);
guidata(gcbo,handles);

% --- Executes on button press in CLGVisSpiral.
function CLGVisSpiral_Callback(hObject, eventdata, handles)
% hObject    handle to CLGVisSpiral (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of CLGVisSpiral
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
% hObject    handle to HSVisSpiral (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of HSVisSpiral
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
% hObject    handle to ZoomChk (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of ZoomChk
axes(handles.MainAxes);
if get(hObject,'Value')
    zoom on
else
    zoom off
end


% --- Executes on button press in CursorChk.
function CursorChk_Callback(hObject, eventdata, handles)
% hObject    handle to CursorChk (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of CursorChk
axes(handles.MainAxes);
if get(hObject,'Value')
    datacursormode on
else
    datacursormode off
end


% --- Executes on button press in Help.
function Help_Callback(hObject, eventdata, handles)
% hObject    handle to Help (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
OFAMMfilePath = mfilename('fullpath');
slashPos = find(OFAMMfilePath == '\');
PathName = OFAMMfilePath(1:slashPos(end));
ReadmePathName = [PathName 'Readme.txt'];
eval(['!notepad ' ReadmePathName])


% --- Executes on button press in AboutOFAMM.
function AboutOFAMM_Callback(hObject, eventdata, handles)
% hObject    handle to AboutOFAMM (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if ~isfield(handles,'AboutFig')
    posMonitor = get(0,'MonitorPositions');
    w = posMonitor(3)/6;
    posAboutFig(3) = w;
    posAboutFig(4) = w;%*posMonitor(3)/posMonitor(4);
    posAboutFig(1) = posMonitor(3)/2-posAboutFig(3)/2;
    posAboutFig(2) = posMonitor(4)/2-0*posAboutFig(4)/2;
    figure('MenuBar','none','ToolBar','none','units',get(0,'units'),'position',posAboutFig);
    
    handles.AboutFig = gcf;
    set(handles.AboutFig,'visible','on','numbertitle','off','Resize','off','name','About OFAMM');
end
AboutOFAMMWin

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
xt = minX + dx*0.1;
yt = minY + dy*0.65;
text(xt,yt,txtStr,'fontsize',9,'FontWeight','bold','fontname','times','HorizontalAlignment','left');

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
