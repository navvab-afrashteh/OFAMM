function handles = NbinSpTraj_fun(handles)

try
    NbinSpTraj_default = 20;
    NbinSpTraj = Str2NumFromHandle(handles.NbinSpTraj,NbinSpTraj_default);
    NbinSpTraj = round(NbinSpTraj);
    set(handles.NbinSpTraj,'string',NbinSpTraj);
catch
    NbinSpTraj = round(handles.NbinSpTraj);
end
handles.TrajHist.NbinSpTraj = NbinSpTraj;