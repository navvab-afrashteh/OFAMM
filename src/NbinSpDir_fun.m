function handles = NbinSpDir_fun(handles)

try
    NbinSpDir_default = 20;
    NbinSpDir = Str2NumFromHandle(handles.NbinSpDir,NbinSpDir_default);
    NbinSpDir = round(NbinSpDir);
    set(handles.NbinSpDir,'string',NbinSpDir);
catch
    NbinSpDir = round(handles.NbinSpDir);
end
handles.HistRose.NbinSpDir = NbinSpDir;
