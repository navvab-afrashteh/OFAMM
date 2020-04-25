function handles = SaveSpeedAngle_fun(handles)


handles = FstartROI_fun(handles);
handles = FendROI_fun(handles);
if isfield(handles,'FstartOFcalculated')
    FstartROI = max([handles.ROI.Fstart, handles.FstartOFcalculated]);
else
    FstartROI = handles.ROI.Fstart;
end
if isfield(handles,'FendOFcalculated')
    FendROI = min([handles.ROI.Fend,handles.FendOFcalculated]);
else
    FendROI = handles.ROI.Fend;
end
try
    set(handles.FstartROI,'string',FstartROI);
    set(handles.FendROI,'string',FendROI);
end
% CLG InstSpAng
if isfield(handles,'InstSpAngcalculated')
    if isfield(handles,'ROI')
        if isfield(handles.ROI,'selected')
            if handles.ROI.selected && isfield(handles,'FstartROIuvCLGcalculated') &&...
                    isfield(handles,'FendROIuvCLGcalculated') && isfield(handles,'ROIuvCLGcalculated')
                if (handles.ROI.Fstart == handles.FstartROIuvCLGcalculated) &&...
                        (handles.ROI.Fend == handles.FendROIuvCLGcalculated) &&...
                        (size(handles.ROI.xy,2) == size(handles.ROIuvCLGcalculated,2))
                    if handles.ROI.xy == handles.ROIuvCLGcalculated
                        saveInstSpAngCLG = 1;
                    else
                        saveInstSpAngCLG = 0; % go and calc inst for CLG
                    end
                else
                    saveInstSpAngCLG = 0; % go and calc inst for CLG
                end
                saveInstSpAngCLG = 1;
            else
                saveInstSpAngCLG = 0; % go and calc inst for CLG
            end
        else
            saveInstSpAngCLG = 0; % go and calc inst for CLG
        end
    end
else
    if isfield(handles,'uvCLGcalculated')
        if  handles.uvCLGcalculated
            %             if ~isfield(handles,'ROI')
            %                 handles = SelectROI(handles.SelectROI, eventdata, handles);
            %             end
            saveInstSpAngCLG = 0; % go and calc inst for CLG
        end
    end
end

% HS InstSpAng
if isfield(handles,'InstSpAngcalculated')
    if isfield(handles,'ROI')
        if isfield(handles.ROI,'selected')
            if handles.ROI.selected && isfield(handles,'FstartROIuvHScalculated') &&...
                    isfield(handles,'FendROIuvHScalculated') && isfield(handles,'ROIuvHScalculated')
                if (handles.ROI.Fstart == handles.FstartROIuvHScalculated) &&...
                        (handles.ROI.Fend == handles.FendROIuvHScalculated) &&...
                        (size(handles.ROI.xy,2) == size(handles.ROIuvCLGcalculated,2))
                    if handles.ROI.xy == handles.ROIuvCLGcalculated
                        saveInstSpAngHS = 1;
                    else
                        saveInstSpAngHS = 0; % go and calc inst for HS
                    end
                else
                    saveInstSpAngHS = 0; % go and calc inst for HS
                end
                saveInstSpAngHS = 1;
            else
                saveInstSpAngHS = 0; % go and calc inst for HS
            end
        else
            saveInstSpAngHS = 0; % go and calc inst for HS
        end
    end
else
    if isfield(handles,'uvHScalculated')
        if  handles.uvHScalculated
            %             if ~isfield(handles,'ROI')
            %                 handles = SelectROI(handles.SelectROI, eventdata, handles);
            %             end
            saveInstSpAngHS = 0; % go and calc inst for HS
        end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CLG TempSpTraj
if isfield(handles,'TempSpTrajLengthcalculated')
    if isfield(handles,'ROI')
        if isfield(handles.ROI,'selected')
            if handles.ROI.selected && isfield(handles,'FstartROITempCLGcalculated') &&...
                    isfield(handles,'FendROITempCLGcalculated') && isfield(handles,'ROITempCLGcalculated')
                if (handles.ROI.Fstart == handles.FstartROITempCLGcalculated) &&...
                        (handles.ROI.Fend == handles.FendROITempCLGcalculated) &&...
                        (size(handles.ROI.xy,2) == size(handles.ROITempCLGcalculated,2))
                    if handles.ROI.xy == handles.ROITempCLGcalculated
                        saveTempSpTrajCLG = 1;
                    else
                        saveTempSpTrajCLG = 0; % go and calc temp for CLG
                    end
                else
                    saveTempSpTrajCLG = 0; % go and calc temp for CLG
                end
                saveTempSpTrajCLG = 1;
            else
                saveTempSpTrajCLG = 0; % go and calc temp for CLG
            end
        else
            saveTempSpTrajCLG = 0; % go and calc temp for CLG
        end
    end
else
    if isfield(handles,'uvCLGcalculated')
        if  handles.uvCLGcalculated
            %             if ~isfield(handles,'ROI')
            %                 handles = SelectROI(handles.SelectROI, eventdata, handles);
            %             end
            saveTempSpTrajCLG = 0; % go and calc temp for CLG
        end
    end
end

% HS TempSpTraj
if isfield(handles,'TempSpTrajLengthcalculated')
    if isfield(handles,'ROI')
        if isfield(handles.ROI,'selected')
            if handles.ROI.selected  && isfield(handles,'FstartROITempHScalculated') &&...
                    isfield(handles,'FendROITempHScalculated') && isfield(handles,'ROITempHScalculated')
                if (handles.ROI.Fstart == handles.FstartROITempHScalculated) &&...
                        (handles.ROI.Fend == handles.FendROITempHScalculated) &&...
                        (size(handles.ROI.xy,2) == size(handles.ROITempHScalculated,2))
                    if handles.ROI.xy == handles.ROITempHScalculated
                        saveTempSpTrajHS = 1;
                    else
                        saveTempSpTrajHS = 0; % go and calc temp for HS
                    end
                else
                    saveTempSpTrajHS = 0; % go and calc temp for HS
                end
                saveTempSpTrajHS = 1;
            else
                saveTempSpTrajHS = 0; % go and calc temp for HS
            end
        else
            saveTempSpTrajHS = 0; % go and calc temp for HS
        end
    end
else
    if isfield(handles,'uvHScalculated')
        if  handles.uvHScalculated
            %             if ~isfield(handles,'ROI')
            %                 handles = SelectROI(handles.SelectROI, eventdata, handles);
            %             end
            saveTempSpTrajHS = 0; % go and calc temp for HS
        end
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% go and calc inst for CLG & HS
if isfield(handles,'uvCLGcalculated') || isfield(handles,'uvHScalculated')
    if  handles.uvCLGcalculated || handles.uvHScalculated
        if saveInstSpAngCLG == 0 || saveInstSpAngHS == 0
            handles = InstSpDir(handles);
            if handles.ROI.selected
                plotflag = 0;
                handles = plotInstSpDir(handles, plotflag);
                handles.InstSpAngcalculated = 1;
            end
        end
    end
end

% go and calc temp for CLG & HS
if isfield(handles,'uvCLGcalculated') || isfield(handles,'uvHScalculated')
    if  handles.uvCLGcalculated || handles.uvHScalculated
        if saveTempSpTrajCLG == 0 || saveTempSpTrajHS == 0
            handles = TempSpTrajLen(handles);
            if handles.ROI.selected
                plotflag = 1;
                handles = plotTempSpTrajLen(handles, plotflag);
                handles.TempSpTrajLengthcalculated = 1;
            end
        end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ~isfield(handles.ROI,'selected')
    handles.ROI.selected = 0;
end
if handles.ROI.selected
    handles = ROISaveNameSpeedAngle_fun(handles);
    FileName_save_SpAngTraj = handles.FileName_save_SpAngTraj;
    
    if isfield(handles,'PathName')
        handles.SavePathName = handles.PathName;
    else
        handles.SavePathName = pwd; % current folder
    end
    % save for CLG
    if isfield(handles,'uvCLGcalculated')
        if  handles.uvCLGcalculated
            if ~exist(handles.SavePathName,'dir')
                mkdir(handles.SavePathName);
            end
            SaveFullFileName = [handles.SavePathName,'\',FileName_save_SpAngTraj,'.mat'];
            if exist(SaveFullFileName)
                load(SaveFullFileName);
            end
            mFileROI_InstSpAng_TempSpTraj = matfile(SaveFullFileName,'Writable',true);
            
            InstSpAng.CLG.uv = handles.InstSpDirCLG;
            InstSpAng.CLG.SpHist.ele = handles.InstSpCLG_ele;
            InstSpAng.CLG.SpHist.bins = handles.InstSpCLG_bins;
            InstSpAng.CLG.AngHist.rho = handles.InstDirCLG_rho;
            InstSpAng.CLG.AngHist.theta = handles.InstDirCLG_theta;
            
            mFileROI_InstSpAng_TempSpTraj.InstSpAng = InstSpAng;
            
            TempSpTraj.CLG.StrLines = handles.TempSpTrajLenCLG;
            TempSpTraj.CLG.TempSpHist.ele = handles.TempSpCLG_ele;
            TempSpTraj.CLG.TempSpHist.bins = handles.TempSpCLG_bins;
            TempSpTraj.CLG.TempLenHist.ele = handles.TempLenCLG_ele;
            TempSpTraj.CLG.TempLenHist.bins = handles.TempLenCLG_bins;
            TempSpTraj.CLG.MaxTempSpHist.ele = handles.MaxTempSpCLG_ele;
            TempSpTraj.CLG.MaxTempSpHist.bins = handles.MaxTempSpCLG_bins;
            
            mFileROI_InstSpAng_TempSpTraj.TempSpTraj = TempSpTraj;
        end
    end
    
    % save for HS
    if isfield(handles,'uvHScalculated')
        if  handles.uvHScalculated
            if ~exist(handles.SavePathName,'dir')
                mkdir(handles.SavePathName);
            end
            SaveFullFileName = [handles.SavePathName,'\',FileName_save_SpAngTraj,'.mat'];
            if exist(SaveFullFileName)
                load(SaveFullFileName);
            end
            mFileROI_InstSpAng_TempSpTraj = matfile(SaveFullFileName,'Writable',true);
            
            InstSpAng.HS.uv = handles.InstSpDirHS;
            InstSpAng.HS.SpHist.ele = handles.InstSpHS_ele;
            InstSpAng.HS.SpHist.bins = handles.InstSpHS_bins;
            InstSpAng.HS.AngHist.rho = handles.InstDirHS_rho;
            InstSpAng.HS.AngHist.theta = handles.InstDirHS_theta;
            
            mFileROI_InstSpAng_TempSpTraj.InstSpAng = InstSpAng;
            
            TempSpTraj.HS.StrLines = handles.TempSpTrajLenHS;
            TempSpTraj.HS.TempSpHist.ele = handles.TempSpHS_ele;
            TempSpTraj.HS.TempSpHist.bins = handles.TempSpHS_bins;
            TempSpTraj.HS.TempLenHist.ele = handles.TempLenHS_ele;
            TempSpTraj.HS.TempLenHist.bins = handles.TempLenHS_bins;
            TempSpTraj.HS.MaxTempSpHist.ele = handles.MaxTempSpHS_ele;
            TempSpTraj.HS.MaxTempSpHist.bins = handles.MaxTempSpHS_bins;
            
            mFileROI_InstSpAng_TempSpTraj.TempSpTraj = TempSpTraj;
        end
    end
    
end

