function handles = TempSpTrajLen(handles)
methods = {'CLG','HS'};
for method = 1:2
    method = methods{method};
    if isfield(handles,['uv',method,'calculated'])
        if eval(['handles.uv',method,'calculated'])
            handles = get_TempSpTrajLen(handles,method);
        end
    end
end