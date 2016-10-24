function handles = calc_save_SourceSink(handles,method)
% source and sink detection
% Detection using Poincare and Jacobian Matrix

Mask = handles.Mask;
FstartOF = handles.FstartOFcalculated;
FendOF = handles.FendOFcalculated;
eval(['uv = handles.uv',method,';']);

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
if ~exist([SSParamsPathName, '\Source-Sink Detection Parameters.txt'],'file')
    generate_source_sink_detection_parameters(SSParamsPathName)
end

fileID = fopen([SSParamsPathName, '\Source-Sink Detection Parameters.txt']);
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
    Contours_Stack, SourceContours_Stack, SinkContours_Stack);

% save the results
% Poincare-Jacobian methods
% "save_source" is a matrix in this form: [x; y; frame]
[d1, d2, d3] = size(SourcePoincareJacobian);
nFrames = d3 +1;
N_all_sources = sum(SourcePoincareJacobian(:));
save_source = zeros(3, N_all_sources);
count = 0;
for idx = 1:nFrames-1
    [r, c] = find(SourcePoincareJacobian(:,:,idx) == 1);
    save_source(:, count+1:count+length(r)) = [c'; r'; idx*ones(1,length(r))+FstartOF-1];
    count = count + length(r);
end

% "save_sink" is a matrix in this form: [x; y; frame]
N_all_sinks = sum(SinkPoincareJacobian(:));
save_sinks = zeros(3, N_all_sinks);
count = 0;
for idx = 1:nFrames-1
    [r, c] = find(SinkPoincareJacobian(:,:,idx) == 1);
    save_sinks(:, count+1:count+length(r)) = [c'; r'; idx*ones(1,length(r))+FstartOF-1];
    count = count + length(r);
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
[X,Y] = meshgrid(1:d2, 1:d1);

count1 = 0;
count2 = 0;

for idx = 1:nFrames-1
    sourceVerified = SourceVerifiedContours_Stack{idx};
    sinkVerified = SinkVerifiedContours_Stack{idx};
    C = Contours_Stack{idx};
    
    for ii = 1:size(sourceVerified,2)
        k = sourceVerified(1,ii);
        n = sourceVerified(2,ii);
        level = sourceVerified(3,ii);
        
        x = C(1,k+1:k+n);
        y = C(2,k+1:k+n);        
        IN = inpolygon(X,Y, x,y);
        area = sum(IN(:));
        
        count1 = count1+1;
        save_contour_source{count1}.xy = [x;y];
        save_contour_source{count1}.frame = idx+FstartOF-1;
        save_contour_source{count1}.strength = level;
        save_contour_source{count1}.size = area;
    end
    
    for ii = 1:size(sinkVerified,2)
        k = sinkVerified(1,ii);
        n = sinkVerified(2,ii);
        level = sinkVerified(3,ii);
        
        x = C(1,k+1:k+n);
        y = C(2,k+1:k+n);
        IN = inpolygon(X,Y, x,y);
        area = sum(IN(:));
        
        count2 = count2+1;
        save_contour_sink{count2}.xy = [x;y];
        save_contour_sink{count2}.frame = idx+FstartOF-1;
        save_contour_sink{count2}.strength = level;
        save_contour_sink{count2}.size = area;
    end
end


eval(['handles.Source',method,'=save_source;']);
eval(['handles.Sink',method,'=save_sinks;']);
eval(['handles.ContourSource',method,'=save_contour_source;']);
eval(['handles.ContourSink',method,'=save_contour_sink;']);

eval(['handles.SS',method,'calculated=1;']); % indicating that source/sink analysis for 'method' has been applieed.

% save source/sink analysis if applicable
saveSS = get(handles.SaveSS,'value');
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
    SaveFullFileName = [handles.SavePathName,'\SourceSinkResults.mat'];
    if exist(SaveFullFileName)
       load(SaveFullFileName);
    end
    eval(['SoSi.',method,'.source = save_source;']);
    eval(['SoSi.',method,'.sink = save_sinks;']);
    eval(['SoSi.',method,'.contour_source = save_contour_source;']);
    eval(['SoSi.',method,'.contour_sink = save_contour_sink;']);
    
    mFileSS = matfile(SaveFullFileName,'Writable',true);
    mFileSS.SoSi = SoSi;
    delete(mFileSS);
end
