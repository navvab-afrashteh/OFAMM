function handles = RunOpticalFlowAnalysisButton(handles)
if ~isfield(handles,'ImgSeqLoaded')
    handles.ImgSeqLoaded = 0;
end
if handles.ImgSeqLoaded
    % CLG
    runCLG = get(handles.runCLG,'value');
    % CLG default values
    CLGparams_default.alpha = 0.03;
    CLGparams_default.ratio = 0.5;
    CLGparams_default.minWidth = round(min(handles.dim1,handles.dim2)*CLGparams_default.ratio/2);
    CLGparams_default.nOuterFPIterations = 7;
    CLGparams_default.nInnerFPIterations = 1;
    CLGparams_default.nSORIterations = 30;
    CLGparams.alpha = Str2NumFromHandle(handles.alphaCLG,CLGparams_default.alpha);
    CLGparams.ratio = Str2NumFromHandle(handles.ratioCLG,CLGparams_default.ratio);
    CLGparams.minWidth = Str2NumFromHandle(handles.minWidthCLG,CLGparams_default.minWidth);
    CLGparams.nOuterFPIterations = Str2NumFromHandle(handles.nOuterFPIterationsCLG,CLGparams_default.nOuterFPIterations);
    CLGparams.nInnerFPIterations = Str2NumFromHandle(handles.nInnerFPIterationsCLG,CLGparams_default.nInnerFPIterations);
    CLGparams.nSORIterations = Str2NumFromHandle(handles.nSORIterationsCLG,CLGparams_default.nSORIterations);
    
    para = [CLGparams.alpha,CLGparams.ratio,CLGparams.minWidth,...
        CLGparams.nOuterFPIterations,CLGparams.nInnerFPIterations,CLGparams.nSORIterations];
    
    % HS
    runHS = get(handles.runHS,'value');
    % HS default values
    HSparams_default.alpha = 0.35;
    HSparams_default.iterations = 2000;
    HSparams.alpha = Str2NumFromHandle(handles.alphaHS,HSparams_default.alpha);
    HSparams.iterations = Str2NumFromHandle(handles.IterationsHS,HSparams_default.iterations);
    
    % TS
    runTS = get(handles.runTS,'value');
    TSparams_default.CorrWinPixels = 3;
    TSparams_default.CorrWinFrames = handles.nFrames;
    TSparams_default.TargetFrames = round(handles.nFrames/2);
    TSparams_default.MaxLag = round(handles.nFrames/20);
    TSparams.CorrWinPixels = Str2NumFromHandle(handles.CorrWinPixelsTS,TSparams_default.CorrWinPixels);
    TSparams.CorrWinFrames = Str2NumFromHandle(handles.CorrWinFramesTS,TSparams_default.CorrWinFrames);
    TSparams.MaxLag = Str2NumFromHandle(handles.MaxLagTS,TSparams_default.MaxLag);
    TSparams.TargetFrames = Str2NumFromHandle(handles.TargetFramesTS,TSparams_default.TargetFrames);
    
    % find first and last frame to calculate OF
    FstartOF_default = 1;
    FendOF_default = handles.nFrames;
    FstartOF = Str2NumFromHandle(handles.FstartOF,FstartOF_default); FstartOF = ceil(FstartOF);
    FendOF = Str2NumFromHandle(handles.FendOF,FendOF_default); FendOF = ceil(FendOF);
    % update fields
    set(handles.FstartOF,'string',num2str(FstartOF));
    set(handles.FendOF,'string',num2str(FendOF));
    
    if runCLG
        handles.uvCLG = zeros(handles.dim1, handles.dim2, FendOF-FstartOF);
    end
    if runHS
        handles.uvHS = zeros(handles.dim1, handles.dim2, FendOF-FstartOF);
    end
    % run OF for CLG or HS if applicable
    if runCLG || runHS
        hWaitBar = waitbar(0,sprintf(' '));
        nFramesTodoOF = FendOF-FstartOF;
        handles.tCLG = 0;
        handles.tHS = 0;
        for idx = FstartOF:FendOF-1
            if runCLG && runHS
                waitbar((idx-FstartOF+1)/nFramesTodoOF,hWaitBar,sprintf('Processing frame %d of %d with CLG and HS',idx-FstartOF+1,nFramesTodoOF));
            elseif runHS
                waitbar((idx-FstartOF+1)/nFramesTodoOF,hWaitBar,sprintf('Processing frame %d of %d with HS',idx-FstartOF+1,nFramesTodoOF));
            elseif runCLG
                waitbar((idx-FstartOF+1)/nFramesTodoOF,hWaitBar,sprintf('Processing frame %d of %d with CLG',idx-FstartOF+1,nFramesTodoOF));
            end
            im1 = handles.ImgSeq(:,:,idx);
            im2 = handles.ImgSeq(:,:,idx+1);
%             [im1, im2] = normalize_two_consequtive_frames(im1,im2,handles.idxMask);
            % CLG
            if runCLG
                tic
                [u, v, ~] = Coarse2FineTwoFrames(im1,im2,para);
                t = toc;
                handles.tCLG = handles.tCLG+t;
                handles.uvCLG(:,:,idx) = (u +1i*v);
            end
            % HS
            if runHS
                tic
                [u, v] = HS(im1,im2,HSparams.alpha,HSparams.iterations);
                t = toc;
                handles.tHS = handles.tHS+t;
                handles.uvHS(:,:,idx) = (u +1i*v);
            end
        end
        handles.FstartOFcalculated = FstartOF;
        handles.FendOFcalculated = FendOF;
        close(hWaitBar);
    end
    
    if runCLG
        handles.uvCLGcalculated = 1;
    end
    if runHS
        handles.uvHScalculated = 1;
    end
    % run TS if applicable
    if runTS
        tic
        handles.uvTS = TS(TSparams, handles);
        handles.tTS = toc;
        handles.uvTScalculated = 1;
    end
    
    % Save section
    if get(handles.SaveOF,'value')
        saveCLG = get(handles.CLGsave,'value');
        saveHS = get(handles.HSsave,'value');
        saveTS = get(handles.TSsave,'value');
        
        % set saving path name
        if saveCLG || saveHS || saveTS
            if isfield(handles,'PathName')
                handles.SavePathName = handles.PathName;
            else
                handles.SavePathName = pwd; % current folder
            end
        end
        
        % save the results
        if runCLG && saveCLG
            if ~exist(handles.SavePathName,'dir')
                mkdir(handles.SavePathName);
            end
            SaveFullFileName = [handles.SavePathName,'\uvResults.mat'];
            mFileuvResults = matfile(SaveFullFileName,'Writable',true);
            mFileuvResults.uvCLG = handles.uvCLG;
            mFileuvResults.FstartOFcalculated =  handles.FstartOFcalculated;
            mFileuvResults.FendOFcalculated = handles.FendOFcalculated;
            mFileuvResults.tCLG = handles.tCLG;
            mFileuvResults.CLGparams = CLGparams;
            delete(mFileuvResults);
            % save CLG parameters
            fileID = fopen([handles.SavePathName, '\CLG Parameters.txt'],'w');
            fprintf(fileID,'%13s\r\n','CLG Parameters');
            fprintf(fileID,'%20s %15s\r\n\r\n','Parameter','Value');
            fprintf(fileID,'%20s %15.3f\r\n','alpha',CLGparams.alpha);
            fprintf(fileID,'%20s %15.2f\r\n','ratio',CLGparams.ratio);
            fprintf(fileID,'%20s %15d\r\n','minWidth',CLGparams.minWidth);
            fprintf(fileID,'%20s %15d\r\n','nOuterFPIterations',CLGparams.nOuterFPIterations);
            fprintf(fileID,'%20s %15d\r\n','nInnerFPIterations',CLGparams.nInnerFPIterations);
            fprintf(fileID,'%20s %15d\r\n','nSORIterations',CLGparams.nSORIterations);
            fclose(fileID);
        end
        if runHS && saveHS
            if ~exist(handles.SavePathName,'dir')
                mkdir(handles.SavePathName);
            end
            SaveFullFileName = [handles.SavePathName,'\uvResults.mat'];
            mFileuvResults = matfile(SaveFullFileName,'Writable',true);
            mFileuvResults.uvHS = handles.uvHS;
            mFileuvResults.FstartOFcalculated =  handles.FstartOFcalculated;
            mFileuvResults.FendOFcalculated = handles.FendOFcalculated;
            mFileuvResults.tHS = handles.tHS;
            mFileuvResults.HSparams = HSparams;
            delete(mFileuvResults);
            % save HS parameters
            fileID = fopen([handles.SavePathName, '\HS Parameters.txt'],'w');
            fprintf(fileID,'%13s\r\n','HS Parameters');
            fprintf(fileID,'%20s %15s\r\n\r\n','Parameter','Value');
            fprintf(fileID,'%20s %15.3f\r\n','alpha',HSparams.alpha);
            fprintf(fileID,'%20s %15d\r\n','iterations',HSparams.iterations);
            fclose(fileID);
        end
        if runTS && saveTS
            if ~exist(handles.SavePathName,'dir')
                mkdir(handles.SavePathName);
            end
            SaveFullFileName = [handles.SavePathName,'\uvResults.mat'];
            mFileuvResults = matfile(SaveFullFileName,'Writable',true);
            mFileuvResults.uvTS = handles.uvTS;
            mFileuvResults.tTS = handles.tTS;
            mFileuvResults.TSparams = TSparams;
            delete(mFileuvResults);
            % save TS parameters
            fileID = fopen([handles.SavePathName, '\TS Parameters.txt'],'w');
            fprintf(fileID,'%13s\r\n','TS Parameters');
            fprintf(fileID,'%20s %15s\r\n\r\n','Parameter','Value');
            fprintf(fileID,'%20s %15d\r\n','CorrWinPixels',TSparams.CorrWinPixels);
            fprintf(fileID,'%20s %15d\r\n','CorrWinFrames',TSparams.CorrWinFrames);
            fprintf(fileID,'%20s %15d\r\n','TargetFrames',TSparams.TargetFrames);
            fprintf(fileID,'%20s %15d\r\n','MaxLag',TSparams.MaxLag);
            fclose(fileID);
        end
    end
end

