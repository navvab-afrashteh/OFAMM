function handles = FendROI_fun(handles)
if isfield(handles,'nFrames')
    FendROI_default = handles.nFrames;
else
    FendROI_default = 1;
end
try
    FendROI = Str2NumFromHandle(handles.FendROI,FendROI_default);
    FendROI = round(FendROI);
    set(handles.FendROI,'string',FendROI);
catch
    FendROI = round(handles.FendROI);
end
handles.ROI.Fend = FendROI;
