function uv = TS(TSparams, handles)

k = TSparams.CorrWinPixels;
corrwin = TSparams.CorrWinFrames;
maxlag = TSparams.MaxLag;
frames = TSparams.TargetFrames;

[dim1, dim2, ~] = size(handles.ImgSeq);

win = CircularWindow(k);

a1 = win(1:(k+1)/2,:);
a1flip = flipud(a1);
[x1, y1] = find(a1flip == 1);

a2 = win(k+1+(k+1)/2:end,:);
[x2, y2] = find(a2 == 1);

a3 = win((k+1)/2+1:-1+k+1+(k+1)/2,1);
[x3, y3] = find(a3 == 1);

a4 = win((k+1)/2+1:-1+k+1+(k+1)/2,end);
[x4, y4] = find(a4 == 1);

points = [k+(k+1)/2+flipud(x2), flipud(y2); (k+1)/2 + flipud(x3), y3; (k+1)/2+1-x1, y1; (k+1)/2 + x4, 2*k+y4];

a=zeros(2*k+1);
for kk=1:length(points)
    a(points(kk,1),points(kk,2))=1;
end

[T1, T2, T3, T4, dir] = Templates(k, points);

uv = zeros(dim1,dim2,length(frames));
p = 1;

hWaitBar = waitbar(0,sprintf(' '));
nFramesTodoOF = length(frames);
frameN = 0;
for frame = frames
    frameN = frameN+1;
    waitbar(frameN/nFramesTodoOF,hWaitBar,sprintf('Processing frame %d with TS',frame));
    for i = k+1:dim1-k-1
        for j = k+1:dim2-k-1
            
            [r2Opt, LagOpt, FR]= OptCorrLagFR(i, j, frame, k, points, corrwin, maxlag, handles);
            Fstar = r2Opt .* LagOpt;
            
            T1star = T1 .* r2Opt;
            T2star = T2 .* r2Opt;
            T3star = T3 .* r2Opt;
            T4star = T4 .* r2Opt;
            
            P1 = (Fstar' * T1star) / (T1star' * T1star);
            if isnan(P1)
                P1 = 0;
            end
            P2 = (Fstar' * T2star) / (T2star' * T2star); 
            if isnan(P2)
                P2 = 0;
            end
            P3 = (Fstar' * T3star) / (T3star' * T3star);
            if isnan(P3)
                P3 = 0;
            end
            P4 = (Fstar' * T4star) / (T4star' * T4star);
            if isnan(P4)
                P4 = 0;
            end
            
            Fplus = P1*T1 + P2*T2 + P3*T3 + P4*T4;
            uv(i,j,p) = (Fplus' * dir);
            
        end
    end
    p = p+1;
end
uv = uv/abs((T1'*dir)^2/(T1'*T1)); % flow in the form of (1/speed)*exp(1i*phi)
uv = uv./(abs(uv).^2); % change the flow to be in the form of (speed)*exp(1i*phi)

% edge pixels
p = 1;
for frame = frames
    for i = 1:dim1
        for j = 1:dim2
            
            if (i<=k || j<=k || i>=dim1-k || j>=dim1-k)
                [r2Opt, LagOpt, FR]= OptCorrLagFR(i, j, frame, k, points, corrwin, maxlag, handles);
                Fstar = r2Opt .* LagOpt;
                
                T1star = T1 .* r2Opt;
                T2star = T2 .* r2Opt;
                T3star = T3 .* r2Opt;
                T4star = T4 .* r2Opt;
                
                P1 = (Fstar' * T1star) / (T1star' * T1star); 
                if isnan(P1)
                    P1 = 0;
                end
                P2 = (Fstar' * T2star) / (T2star' * T2star); 
                if isnan(P2)
                    P2 = 0;
                end
                P3 = (Fstar' * T3star) / (T3star' * T3star);
                if isnan(P3)
                    P3 = 0;
                end
                P4 = (Fstar' * T4star) / (T4star' * T4star);
                if isnan(P4)
                    P4 = 0;
                end
                
                Fplus = P1*T1 + P2*T2 + P3*T3 + P4*T4;                
                uv(i,j,p) = (Fplus' * (dir .* (r2Opt~=0)));
                uv(i,j,p) = uv(i,j,p)/abs((T1'*(dir .* (r2Opt~=0)))^2/(T1'*(T1 .* (r2Opt~=0))));
                uv(i,j,p) = uv(i,j,p)/(abs(uv(i,j,p))^2);
            end
        end
    end
    p = p+1;
end
close(hWaitBar);


