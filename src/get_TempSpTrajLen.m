function handles = get_TempSpTrajLen(handles,method)
eventdata = [];
if isfield(handles,['uv',method,'calculated'])
    if eval(['handles.uv',method,'calculated'])
%         [~, flag] = getimage(handles.MainAxes); % check if there is an image in the main axes
% %         if flag
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
                sx = xy(1,:);
                sy = xy(2,:);
                startFrame = FstartROI;
                endFrame = FendROI;
                offsetFrame = 1;
                if startFrame < endFrame
                    output = Velocity_Profile(eval(['handles.uv',method]),startFrame,endFrame,offsetFrame,sx,sy);
                    
                    eval(['handles.TempSpTrajLen',method,' = output;']);
                    eval(['handles.FstartROITemp',method,'calculated = FstartROI;']);
                    eval(['handles.FendROITemp',method,'calculated = FendROI;']);
                    eval(['handles.ROITemp',method,'calculated = xy;']);
                end
            end
%         end
    end
end
