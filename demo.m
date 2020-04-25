clc
clear
close all

%% loadin a .tiff image sequence. 
% change these for your data
handles.PathName = 'E:\Github Repositories\OFAMM-Sample-Data\Real_WFI_Data\AuditoryStim_VSDI\';
FullFileName = sprintf('%s//%s',handles.PathName,'ImgSeq.tif');
handles.ImgSeq = imreadalltiff(FullFileName);
[handles.dim1,handles.dim2,handles.nFrames] = size(handles.ImgSeq);
handles.ImgSeqLoaded = 1;
handles.uvCLGcalculated = 0;
handles.uvHScalculated = 0;
handles.uvTScalculated = 0;
% to use a .raw data use ImgSeq = imreadallraw(FullFileName,xdim,ydim,nFrames,'*float32');
% to use a .raw data use ImgSeq = load(FullFileName);

%% load the craniotomy window mask
MaskFileName = sprintf('%s//%s',handles.PathName,'Mask.tif');
handles.Mask = imread(MaskFileName)>0;
[handles.rMask,handles.cMask] = find(handles.Mask);
handles.idxMask = sub2ind(size(handles.Mask),handles.rMask,handles.cMask);
% I create an ROI as well. You can modify or load your desired ROI(s).
% Note that ROI should be logical.
handles.ROI.BW = false(size(handles.Mask)); handles.ROI.BW(80:100,30:40) = true;
[r,c] = find(handles.ROI.BW);
handles.ROI.xy = [c';r'];
handles.ROI.selected = 1;

%% Plot the intensity vs frames for the average values in the ROI
handles.ROI.Fstart = 1;
handles.ROI.Fend = handles.nFrames; % you can change these as desired
FstartROI = handles.ROI.Fstart;
FendROI = handles.ROI.Fend;
handles.IntensitySig = zeros(1,FendROI-FstartROI+1);
for idx = FstartROI:FendROI
    thisframe = handles.ImgSeq(:,:,idx);
    handles.IntensitySig(idx-FstartROI+1) = mean(thisframe(handles.ROI.BW));
end
figure(101);
handles.IntensityFig = gcf;
set(handles.IntensityFig,'visible','on','numbertitle','off','name','Intensity vs. Time for selected ROI');
set(handles.IntensityFig,'units','norm','pos',[0.01,0.6,.3,.3])
plot(FstartROI:FendROI, handles.IntensitySig)
xlabel('frame number')
ylabel('Mean Intensity')

%% optical flow analysis
% run for HS and CLG then save
handles.SaveOF = 1;
% CLG params
handles.runCLG = 1;
handles.saveCLG = 1;
handles.CLGparams.alpha = 0.03;
handles.CLGparams.ratio = 0.5;
handles.CLGparams.minWidth = round(min(handles.dim1,handles.dim2)*handles.CLGparams.ratio/2);
handles.CLGparams.nOuterFPIterations = 7;
handles.CLGparams.nInnerFPIterations = 1;
handles.CLGparams.nSORIterations = 30;
% HS params
handles.runHS = 1;
handles.saveHS = 1;
handles.HSparams.alpha = 0.35;
handles.HSparams.iterations = 2000;
% specify frames that you want to run optical flow on
handles.FstartOF = 1;
handles.FendOF = handles.nFrames;
% run optical flow and save 
handles = RunOpticalFlowAnalysisButton(handles);

%% Source-Sink analysis
handles.SaveSS = 1;
method = 'CLG';
handles = calc_save_SourceSink(handles,method);
method = 'HS';
handles = calc_save_SourceSink(handles,method);

%{
The results are as follows: handles.SourceCLG, handles.SinkCLG,
handles.SourceHS, and handles.SinkHS show the frame number, center of point
x and y. [framNum; x; y]. The size of each of them depends on how many
source or sinks with that method were identified. size = 3 x Num source or sink
3 is the length of [framNum; x; y]

handles.ContourSourceCLG and other similar params have the info about the
source or sink and their contour x and y. Length of
handles.ContourSourceCLG is the same as size(handles.SourceCLG,2). For each
source there is a cell in handles.ContourSourceCLG and it contains xy of
contour, frame number, strength of the source and the pixel size of it.
%}

%% Generate stats on the instantaneous speed and angle and the trajectories 
% and its temporal speed and trajectory length. These are for the selected
% ROI. ROI is selected on line 24.
plotflag = 1;
handles.FstartROI = 40;
handles.FendROI = 90;
% instantaneous speed and angle
handles.NbinSpDir = 15;
handles = InstSpDir(handles);
handles = plotInstSpDir(handles, plotflag);
handles.InstSpAngcalculated = 1;
% Traj and temporal speed
handles.NbinSpTraj = 15;
handles = TempSpTrajLen(handles);
handles = plotTempSpTrajLen(handles, plotflag);
handles.TempSpTrajLengthcalculated = 1;
% saving
handles.ROISaveNameSpeedAngle = 'ROI1'; % name your ROI to save
handles = SaveSpeedAngle_fun(handles);
%{
In the folder there will be a file called ROI1.mat which contains all the
info of instantaneous speed and angle for the pixels in that ROI and the
bar graph distribution of speed and the rose plot of angles.

InstSpAng.CLG.uv : complex velocity values for all the pixels and selected
frames. size: Num pixels in ROI x Num selected frames.
InstSpAng.CLG.SpHist : Speed histogram info. Contains ele (histogram values)
and bins (bins' centers)
InstSpAng.CLG.AngHist : angle rose plot info. contains rho (length for each
theta) and theta (center of each bin)

It alos has info on the trajectory and temporal speed for each pixel in the
ROI. In each of the followings there is a cell for each pixel.

TempSpTraj.CLG.MaxTempSpHist : maximum temporal speed histogram info. contains ele and bins.
TempSpTraj.CLG.TempLenHist : trajectory length historgram info. contains ele and bins.
TempSpTraj.CLG.TempSpHist : average temporal speed hist info. contains ele and bins.

TempSpTraj.CLG.StrLines : all the info about the streamlines (trajectory)
for each pixel.

TempSpTraj.CLG.StrLines.tspeed: time stamps for the trajectory of each
pixel. the info for each pixel is in a cell. These are basically integers
from first to last desired frames.

TempSpTraj.CLG.StrLines.speed:  speed values for the trajectory of each
pixel. the info for each pixel is in a cell. These are basically speeds
from first to last desired frames.

TempSpTraj.CLG.StrLines.str: trajectory location [x, y, frameNum] for each
pixel. Info for each pixel is in a cell. for each pixel there is a matrix
of the size Num frames x 3.
%}