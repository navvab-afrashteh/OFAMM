function handles = InstSpDir(handles)
%ViewInstSpeedAngle
methods = {'CLG','HS'};
for method = 1:2
    method = methods{method};
    if isfield(handles,['uv',method,'calculated'])
        if eval(['handles.uv',method,'calculated'])
            handles = get_InstSpDir(handles,method);
        end
    end
end