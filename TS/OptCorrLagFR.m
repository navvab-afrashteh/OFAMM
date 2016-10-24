function [r2Opt, LagOpt, FR]= OptCorrLagFR(i, j, frame, k, points, corrwin, maxlag, handles)

% Finding the maximum value of correlation between pixel (i,j) and each 
% pixel in the window 'win' surrondimg it, by finding the optimal value for
% time shift

[d1,d2,nFrames] = size(handles.ImgSeq);

if corrwin >= nFrames
    startFrame = 1;
    stopFrame = nFrames;
else
    w = 2*round((corrwin-1)/2);
    startFrame = max(1,frame - w/2);
    stopFrame = min(nFrames,frame + w/2);
    if startFrame == 1
        stopFrame = w;
    end
    if stopFrame == nFrames
        startFrame = nFrames - w + 1;
    end
end

row = points(:,1);
col = points(:,2);
x = row+i-k-1;
y = col+j-k-1;
len = length(points);

ij = handles.ImgSeq(i,j,startFrame:stopFrame);
ij = ij(:);

r2Opt = zeros(2*len,1);
LagOpt = zeros(2*len,1);
for idx = 1:len
    if x(idx)>0 && y(idx)>0 && x(idx)<=d1 && y(idx)<=d2 
        xy = handles.ImgSeq(x(idx) ,y(idx), startFrame:stopFrame);
        [r, lags] = xcorr(ij, xy, maxlag, 'coef');
        [r2Opt(idx,1), lagIndex] = max(r);
        LagOpt(idx,1) = -lags(lagIndex);
    else
        r2Opt(idx,1) = 0;
        LagOpt(idx,1) = 0;
    end
end

xc = circshift(x,[-1, 0]);
yc = circshift(y,[-1, 0]);
for idx = len+1:2*len
    if x(idx-len)>0 && y(idx-len)>0 && xc(idx-len)>0 && yc(idx-len)>0 &&...
            x(idx-len)<=d1 && y(idx-len)<=d2 && xc(idx-len)<=d1 && yc(idx-len)<=d2 
        
        ij = handles.ImgSeq(x(idx-len), y(idx-len), startFrame:stopFrame);
        ij = ij(:);
        xy = handles.ImgSeq(xc(idx-len),yc(idx-len),startFrame:stopFrame);
        xy = xy(:);
        [r, lags] = xcorr(ij, xy, maxlag, 'coef');
        [r2Opt(idx,1), lagIndex] = max(r);
        LagOpt(idx,1) = -lags(lagIndex);
    else
        r2Opt(idx,1) = 0;
        LagOpt(idx,1) = 0;
    end
end
FR = sum(r2Opt)/len; % Flow Reliability
