function handles = FstartROI_fun(handles)
FstartROI_default = 1;
try
    FstartROI = Str2NumFromHandle(handles.FstartROI,FstartROI_default);
    FstartROI = round(FstartROI);
    set(handles.FstartROI,'string',FstartROI);
catch
    FstartROI = round(handles.FstartROI);
end
handles.ROI.Fstart = FstartROI;
