function handles = calc_save_SourceSink(handles,method)
% source and sink detection
% Detection using Poincare and Jacobian Matrix

Mask = handles.Mask;
FstartOF = handles.FstartOFcalculated;
FendOF = handles.FendOFcalculated;
eval(['uv = handles.uv',method,';']);
uv = uv(:,:,FstartOF:FendOF-1);
% source and sink detection
% Detection using Poincare and Jacobian Matrix
[SourcePoincareJacobian, SinkPoincareJacobian] = SourceSink_PoincareJacobian(uv, Mask, method);

%Dection using gradient of vector field
% Source-Sink Detection Parameters
if isfield(handles,'PathName')
    SSParamsPathName = handles.PathName;
else
    SSParamsPathName = pwd;
end
if ~exist(SSParamsPathName,'dir')
    mkdir(SSParamsPathName);
end
if ~exist([SSParamsPathName, '\Source-Sink-Spiral Detection Parameters.txt'],'file')
    generate_source_sink_detection_parameters(SSParamsPathName)
end

fileID = fopen([SSParamsPathName, '\Source-Sink-Spiral Detection Parameters.txt']);
textscan(fileID,'%s %s %s',1);
textscan(fileID,'%s %s',1);
parameters = textscan(fileID,'%s %f',7);
for k = 1:length(parameters{1})
    eval([parameters{1,1}{k}, '=', 'parameters{1,2}(k);']);
end
fclose(fileID);

[SourceContours_Stack, SinkContours_Stack, Contours_Stack] = ...
    SourceSink_Contours(uv, Mask, Nmin, Lmin_source, Lmax_sink);

[SourcePoincareJacobian, SinkPoincareJacobian,...
    SourceVerifiedContours_Stack, SinkVerifiedContours_Stack] = ...
    Confirm_Contours_FPs(SourcePoincareJacobian, SinkPoincareJacobian, ...
    Contours_Stack, SourceContours_Stack, SinkContours_Stack,...
    Nnested_source, Nnested_sink);

% save the results
% Poincare-Jacobian methods
% "save_source" is a matrix in this form: [x; y; frame]
[d1, d2, d3] = size(SourcePoincareJacobian);
nFrames = d3 +1;
N_all_sources = sum(SourcePoincareJacobian(:)==2);
N_all_sources_spiral = sum(SourcePoincareJacobian(:)==1);
save_source = zeros(3, N_all_sources);
save_source_spiral = zeros(3, N_all_sources_spiral);
countN = 0;
countS = 0;
for idx = 1:nFrames-1
    % spiral sources
    [r, c] = find(SourcePoincareJacobian(:,:,idx) == 1);
    save_source_spiral(:, countS+1:countS+length(r)) = [c'; r'; idx*ones(1,length(r))+FstartOF-1];
    countS = countS + length(r);
    % node sources
    [r, c] = find(SourcePoincareJacobian(:,:,idx) == 2);
    save_source(:, countN+1:countN+length(r)) = [c'; r'; idx*ones(1,length(r))+FstartOF-1];
    countN = countN + length(r);
end

% "save_sink" is a matrix in this form: [x; y; frame]
N_all_sinks = sum(SinkPoincareJacobian(:)==2);
N_all_sinks_spiral = sum(SinkPoincareJacobian(:)==1);
save_sinks = zeros(3, N_all_sinks);
save_sinks_spiral = zeros(3, N_all_sinks_spiral);
countN = 0;
countS = 0;
for idx = 1:nFrames-1
    % spiral sources
    [r, c] = find(SinkPoincareJacobian(:,:,idx) == 1);
    save_sinks_spiral(:, countS+1:countS+length(r)) = [c'; r'; idx*ones(1,length(r))+FstartOF-1];
    countS = countS + length(r);
    % node sources
    [r, c] = find(SinkPoincareJacobian(:,:,idx) == 2);
    save_sinks(:, countN+1:countN+length(r)) = [c'; r'; idx*ones(1,length(r))+FstartOF-1];
    countN = countN + length(r);
end

% Divergence method (contours)
% "save_contour_source" and "save_contour_sink" are cells with these fields:
% xy, frame, strength, and area

N_all_contour_sources = 0;
N_all_contour_sinks = 0;
for idx = 1:nFrames-1
    sourceVerified = SourceVerifiedContours_Stack{idx};
    sinkVerified = SinkVerifiedContours_Stack{idx};
    
    N_all_contour_sources = N_all_contour_sources + size(sourceVerified,2);
    N_all_contour_sinks = N_all_contour_sinks + size(sinkVerified,2);
end
save_contour_source = cell(1,N_all_contour_sources);
save_contour_sink = cell(1,N_all_contour_sinks);

save_contour_source = cell(1,N_all_sources);
save_contour_source_spiral = cell(1,N_all_sources_spiral);
save_contour_sink = cell(1,N_all_sinks);
save_contour_sink_spiral = cell(1,N_all_sinks_spiral);

[X,Y] = meshgrid(1:d2, 1:d1);

count1 = 0; count1sp = 0;
count2 = 0; count2sp = 0;

for idx = 1:nFrames-1
    sourceVerified = SourceVerifiedContours_Stack{idx};
    sinkVerified = SinkVerifiedContours_Stack{idx};
    C = Contours_Stack{idx};
    
    S = SourcePoincareJacobian(:,:,idx);
    S = S(S>0);
    for ii = 1:size(sourceVerified,2)
        k = sourceVerified(1,ii);
        n = sourceVerified(2,ii);
        level = sourceVerified(3,ii);
        
        x = C(1,k+1:k+n);
        y = C(2,k+1:k+n);        
        IN = inpolygon(X,Y, x,y);
        area = sum(IN(:));
        
        if S(ii) == 2
            count1 = count1+1;
            save_contour_source{count1}.xy = [x;y];
            save_contour_source{count1}.frame = idx+FstartOF-1;
            save_contour_source{count1}.strength = level;
            save_contour_source{count1}.size = area;
        end
        if S(ii) == 1
            count1sp = count1sp+1;
            save_contour_source_spiral{count1sp}.xy = [x;y];
            save_contour_source_spiral{count1sp}.frame = idx+FstartOF-1;
            save_contour_source_spiral{count1sp}.strength = level;
            save_contour_source_spiral{count1sp}.size = area;
        end
    end
    
    S = SinkPoincareJacobian(:,:,idx);
    S = S(S>0);
    for ii = 1:size(sinkVerified,2)
        k = sinkVerified(1,ii);
        n = sinkVerified(2,ii);
        level = sinkVerified(3,ii);
        
        x = C(1,k+1:k+n);
        y = C(2,k+1:k+n);
        IN = inpolygon(X,Y, x,y);
        area = sum(IN(:));
        
        if S(ii) == 2
            count2 = count2+1;
            save_contour_sink{count2}.xy = [x;y];
            save_contour_sink{count2}.frame = idx+FstartOF-1;
            save_contour_sink{count2}.strength = level;
            save_contour_sink{count2}.size = area;
        end
        if S(ii) == 1
            count2sp = count2sp+1;
            save_contour_sink_spiral{count2sp}.xy = [x;y];
            save_contour_sink_spiral{count2sp}.frame = idx+FstartOF-1;
            save_contour_sink_spiral{count2sp}.strength = level;
            save_contour_sink_spiral{count2sp}.size = area;
        end
    end
end
eval(['handles.Source',method,'=save_source;']);
eval(['handles.Sink',method,'=save_sinks;']);
eval(['handles.ContourSource',method,'=save_contour_source;']);
eval(['handles.ContourSink',method,'=save_contour_sink;']);
eval(['handles.SourceSpiral',method,'=save_source_spiral;']);
eval(['handles.SinkSpiral',method,'=save_sinks_spiral;']);
eval(['handles.ContourSourceSpiral',method,'=save_contour_source_spiral;']);
eval(['handles.ContourSinkSpiral',method,'=save_contour_sink_spiral;']);
eval(['handles.SS',method,'calculated=1;']); % indicating that source/sink analysis for 'method' has been applieed.

% save source/sink analysis if applicable
SaveSS = 0;
try saveSS = get(handles.SaveSS,'value'); catch; try saveSS = handles.SaveSS; end; end
if saveSS
    % set saving path name
    if isfield(handles,'PathName')
        handles.SavePathName = handles.PathName;
    else
        handles.SavePathName = pwd; % current folder
    end
    % save the results
    if ~exist(handles.SavePathName,'dir')
        mkdir(handles.SavePathName);
    end
    SaveFullFileName = [handles.SavePathName,'\SourceSinkSpiralResults.mat'];
    if exist(SaveFullFileName)
       load(SaveFullFileName);
    end
    eval(['SoSi.',method,'.source = save_source;']);
    eval(['SoSi.',method,'.sink = save_sinks;']);
    eval(['SoSi.',method,'.contour_source = save_contour_source;']);
    eval(['SoSi.',method,'.contour_sink = save_contour_sink;']);
    
    eval(['SoSi.',method,'.source_spiral = save_source_spiral;']);
    eval(['SoSi.',method,'.sink_spiral = save_sinks_spiral;']);
    eval(['SoSi.',method,'.contour_source_spiral = save_contour_source_spiral;']);
    eval(['SoSi.',method,'.contour_sink_spiral = save_contour_sink_spiral;']);
    
    mFileSS = matfile(SaveFullFileName,'Writable',true);
    mFileSS.SoSi = SoSi;
    delete(mFileSS);
end
