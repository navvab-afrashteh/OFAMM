function handles = get_InstSpDir(handles,method)
eventdata=[];
if isfield(handles,['uv',method,'calculated'])
    if eval(['handles.uv',method,'calculated'])
        
            if isfield(handles,'ROI')
                if isfield(handles.ROI,'selected')
                    if ~handles.ROI.selected
                        [~, flag] = getimage(handles.MainAxes); % check if there is an image in the main axes
                        if flag
                            handles = SelectROI(handles.SelectROI, eventdata, handles);
                        end
                    end
                else
                    [~, flag] = getimage(handles.MainAxes); % check if there is an image in the main axes
                    if flag
                        handles = SelectROI(handles.SelectROI, eventdata, handles);
                    end
                end
            else
                handles = SelectROI(handles.SelectROI, eventdata, handles);
            end
            
            if handles.ROI.selected
                handles = FstartROI_fun(handles);
                handles = FendROI_fun(handles);
                FstartROI = max([handles.ROI.Fstart, handles.FstartOFcalculated]);
                FendROI = min([handles.ROI.Fend,handles.FendOFcalculated]);
                try
                    set(handles.FstartROI,'string',FstartROI);
                    set(handles.FendROI,'string',FendROI);
                end
                xy = handles.ROI.xy;
                xyind = sub2ind([handles.dim1, handles.dim2],xy(2,:),xy(1,:));
                nPoints = size(xy,2);
                eval(['handles.InstSpDir',method,' = zeros(nPoints,FendROI-FstartROI);']);
                for idx = FstartROI:FendROI-1
                    eval(['thisframe = handles.uv',method,'(:,:,idx);']);
                    eval(['handles.InstSpDir',method,'(:,idx-FstartROI+1) = thisframe(xyind);']);
                end
                eval(['handles.FstartROIuv',method,'calculated = FstartROI;']);
                eval(['handles.FendROIuv',method,'calculated = FendROI;']);
                eval(['handles.ROIuv',method,'calculated = xy;']);
            end
%         end
    end
end
