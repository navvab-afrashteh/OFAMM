function handles = ROISaveNameSpeedAngle_fun(handles)
try
    handles.FileName_save_SpAngTraj = get(handles.ROISaveNameSpeedAngle,'string');
catch
    handles.FileName_save_SpAngTraj = handles.ROISaveNameSpeedAngle;
end
try
    handles.FileName_save_SpAngTraj = strtrim(handles.FileName_save_SpAngTraj);
catch
    handles.FileName_save_SpAngTraj = '';
end
if strcmp(handles.FileName_save_SpAngTraj,'')
    handles.FileName_save_SpAngTraj = 'ROI_InstSpAng_TempSpTraj';
    try
        set(handles.ROISaveNameSpeedAngle,'string',handles.FileName_save_SpAngTraj);
    end
end
