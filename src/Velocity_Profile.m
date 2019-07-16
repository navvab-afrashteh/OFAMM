function output = Velocity_Profile (uv,startFrame,endFrame,offsetFrame,sx,sy,varargin)

% uv = uv(:,:,startFrame+offsetFrame-1:endFrame+offsetFrame-2);
[d1, d2, d3] = size(uv);
uv = cat(3, zeros(d1,d2),uv);

idxS = startFrame;
idxE = endFrame;
st = idxS*ones(size(sx));
nFrames = idxE-idxS+1;

U = real(uv);
V = imag(uv);
W = ones(size(U));

[X,Y,T] = meshgrid(1:size(U,1), 1:size(U,2), 1:size(U,3));
XYT = stream3(X,Y,T,U,V,W,sx,sy,st);

%% calc estimated velocity 
ti = cell(1,length(XYT));
est_velocity_of_propagation = cell(1,length(XYT));
for point = 1:length(XYT)
    xyt_OF = XYT{point};
    x_OF_all = xyt_OF(:, 1);
    y_OF_all = xyt_OF(:, 2);
    t_OF_all = xyt_OF(:, 3);
    
    x_OF = zeros(nFrames,1);
    y_OF = zeros(nFrames,1);
    t_OF = zeros(nFrames,1);
    p=1;
    for idx = idxS:idxE
        [~, ii] = min(abs(t_OF_all - idx));
        x_OF(p) = x_OF_all(ii);
        y_OF(p) = y_OF_all(ii);
        t_OF(p) = t_OF_all(ii)+offsetFrame-1;
        p=p+1;
    end
    XYT{point}  = [x_OF, y_OF, t_OF];
    
    % calculate velocity of propagation
    est_dist_of_propagation = sqrt(diff(x_OF).^2 + diff(y_OF).^2);
    est_time_of_propagation = diff(t_OF);
    est_velocity_of_propagation{point} = (est_dist_of_propagation./est_time_of_propagation);
    ti{point} = t_OF(1:end-1);
end

output.tspeed = ti;
output.speed = est_velocity_of_propagation;
output.str = XYT;


