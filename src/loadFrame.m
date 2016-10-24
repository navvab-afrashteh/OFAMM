
function loadFrame(handles,frameNumber)
% get min and max disp values
handles = MinDispVal_Callback(handles.MinDispVal, eventdata, handles);
handles = MaxDispVal_Callback(handles.MaxDispVal, eventdata, handles);
MinVal = handles.MinVal;
MaxVal = handles.MaxVal;
% imagesc frame and mask
if isfield(handles,'ImgSeq')
    if handles.ImgSeqLoaded
        if frameNumber > 0 && frameNumber <= handles.nFrames
            if get(handles.MaskChk,'value')
                thisframe = handles.ImgSeq(:,:,frameNumber).*handles.Mask;
            else
                thisframe = handles.ImgSeq(:,:,frameNumber);
            end
            cla(handles.MainAxes);
            imagesc(thisframe,[MinVal,MaxVal],'parent',handles.MainAxes);
            colormap jet;
            hold on;
            set(handles.FrameNumber,'String',sprintf('%d',frameNumber));
        end
    end
end
% show optical flow
if get(handles.VectorFieldVisChk,'value')
    if get(handles.CLGVisOF,'value'); methodOF = 'CLG'; end
    if get(handles.HSVisOF,'value'); methodOF = 'HS'; end
    
    if eval(['isfield(handles,','''uv',methodOF,'calculated''))'])
        if eval(['handles.uv',methodOF,'calculated'])
            if frameNumber > 0 && frameNumber <= handles.nFrames &&...
                    frameNumber >= handles.FstartOFcalculated &&...
                    frameNumber < handles.FendOFcalculated
                
                frameNumberOF = frameNumber-handles.FstartOFcalculated+1;
                eval(['thisuv = handles.uv',methodOF,'(:,:,frameNumberOF);']);
                u = real(thisuv(:,:,frameNumber));
                v = imag(thisuv(:,:,frameNumber));
                [gridXDown,gridYDown,OpticalFlowDown_u] = mean_downsample(u,3);
                [~,~,OpticalFlowDown_v]                 = mean_downsample(v,3);
                set(handles.figure1,'currentAxes',handles.MainAxes);
                quiver(gridXDown,gridYDown,OpticalFlowDown_u(:,:,1),OpticalFlowDown_v(:,:,1),2.5,'color','k','LineWidth',0.75);
                hold on;
            end
        end
    end
end
% display source/sink
markersize = 10;
linewidth = 2;
if get(handles.SSChk,'value')
    if get(handles.CLGVisSS,'value'); methodSS = 'CLG'; end
    if get(handles.HSVisSS,'value'); methodSS = 'HS'; end
    
    if eval(['isfield(handles,','''SS',methodSS,'calculated''))'])
        if eval(['handles.SS',methodSS,'calculated'])
            eval(['source = handles.Source',method,';']);
            eval(['sinks = handles.Sink',method,';']);
            eval(['contour_source = handles.ContourSource',method,';']);
            eval(['contour_sink = handles.ContourSink',method,';']);
            %find sources and display
            [~, sourceidx] = find(source(3,:) == frameNumber);
            for idx = sourceidx
                set(handles.figure1,'currentAxes',handles.MainAxes);
                plot(source(1,idx),source(2,idx),'r*','markersize',markersize)
                hold on;
                plot(contour_source{idx}.xy(1,:),contour_source{idx}.xy(2,:),'r.-','linewidth',linewidth)
                hold on
            end
            %find sinks and display
            [~, sinkidx] = find(sink(3,:) == frameNumber);
            for idx = sinkidx
                set(handles.figure1,'currentAxes',handles.MainAxes);
                plot(sink(1,idx),sink(2,idx),'wo','markersize',markersize)
                hold on;
                plot(contour_sink{idx}.xy(1,:),contour_sink{idx}.xy(2,:),'w.-','linewidth',linewidth)
                hold on
            end
        end
    end
end

% display trajectories
if get(handles.Trajectory,'value')
    handles.Trajectory.currVal = 1;
    if get(handles.CLGVisTraj,'value'); methodTraj = 'CLG'; end
    if get(handles.HSVisTraj,'value'); methodTraj = 'HS'; end
    
    if eval(['isfield(handles,','''uv',methodTraj,'calculated''))'])
        if eval(['handles.uv',methodTraj,'calculated'])
            if ~isfield(handles,'Trajectory.preVal')
                handles.Trajectory.preVal = 0;
            end
            
            if handles.Trajectory.preVal
                FstartROI = handles.Trajectory.Fstart;
                FstartROI = max([FstartROI, handles.FstartOFcalculated]);
                FendROI = min([frameNumber,handles.FendOFcalculated]);
                
                xy = handles.ROI.xy;
                sx = xy(1,:);
                sy = xy(2,:);
                startFrame = FstartROI - handles.FstartOFcalculated +1;
                endFrame = FendROI - handles.FstartOFcalculated +1;
                offsetFrame = handles.FstartOFcalculated;
                if startFrame < endFrame
                    output = Velocity_Profile(eval(['handles.uv',methodTraj]),startFrame,endFrame,offsetFrame,sx,sy);
                    % display
                    set(handles.figure1,'currentAxes',handles.MainAxes);
                    str = streamline(output.str);
                    hold on;
                    handles.Trajectory.preVal = 1;
                    handles.Trajectory.Fstart = FstartROI;
                end
            else
                if isfield(handles,'ROI')
                    if isfield(handles.ROI,'selected')
                        if ~handles.ROI.selected
                            handles = SelectROI(handles.SelectROI, eventdata, handles);
                        end
                    end
                end
                
                if handles.ROI.selected
                    handles = FstartROI_Callback(handles.FstartROI, eventdata, handles);
                    FstartROI = max([handles.ROI.Fstart, handles.FstartOFcalculated]);
                    FendROI = min([frameNumber,handles.FendOFcalculated]);
                    set(handles.FstartROI,'string',FstartROI);
                    
                    xy = handles.ROI.xy;
                    sx = xy(1,:);
                    sy = xy(2,:);
                    startFrame = FstartROI - handles.FstartOFcalculated +1;
                    endFrame = FendROI - handles.FstartOFcalculated +1;
                    offsetFrame = handles.FstartOFcalculated;
                    if startFrame < endFrame
                        output = Velocity_Profile(eval(['handles.uv',methodTraj]),startFrame,endFrame,offsetFrame,sx,sy);
                        % display
                        set(handles.figure1,'currentAxes',handles.MainAxes);
                        str = streamline(output.str);
                        hold on
                        handles.Trajectory.preVal = 1;
                        handles.Trajectory.Fstart = FstartROI;
                    end
                end
            end
        end
    end
end

