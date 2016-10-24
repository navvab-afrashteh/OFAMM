function [SourceContours_Stack, SinkContours_Stack, Contours_Stack] = ...
    SourceSink_Contours(uv, Mask, Nmin, Lmin_source, Lmax_sink)

% Nmin: minimum number of points as the contour size
% Lmin_source: minimum source level
% Lmax_sink: maximum sink level

% Finding Contours and Refining them
nFrames = size(uv,3)+1;

SourceContours_Stack = cell(1,nFrames-1);
SinkContours_Stack = cell(1,nFrames-1);
Contours_Stack = cell(1,nFrames-1);

for idx = 1:nFrames-1
    u = real(uv(:,:,idx));
    v = imag(uv(:,:,idx));
    
    [source, sink, C] = refiningContours(u,v,Mask,Nmin, Lmin_source, Lmax_sink);
    
    SourceContours_Stack{idx} = source;
    SinkContours_Stack{idx} = sink;
    Contours_Stack{idx} = C;
end


