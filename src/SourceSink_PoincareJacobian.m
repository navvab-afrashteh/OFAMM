function [SourcePoincareJacobian, SinkPoincareJacobian] = SourceSink_PoincareJacobian(uv, Mask, method)
% this is July 19 2016 to check waitbar
% source and sink detection
% Poincare Fixed-Points Categorization

PoincareSinkSource = zeros(size(uv));
PoincareSaddle = zeros(size(uv));
FixedPointsPoincare = zeros(size(uv));
hWaitBar = waitbar(0,sprintf(' '));
pos = get(hWaitBar,'position');
pos(2) = pos(2)+pos(4)*1.75;
set(hWaitBar,'position',pos);

for idx = 1:size(uv,3)
    waitbar(idx/size(uv,3),hWaitBar,sprintf(['PoincareJacobian on ', method,', frames ... %d of %d ... '],idx,size(uv,3)));
    [PoincareSinkSource(:,:,idx), PoincareSaddle(:,:,idx)] = poincare_index(uv(:,:,idx));
    PoincareSinkSource(:,:,idx) = PoincareSinkSource(:,:,idx) .* Mask;
    PoincareSaddle(:,:,idx) = PoincareSaddle(:,:,idx) .* Mask;
    
    FixedPointsPoincare(:,:,idx) = PoincareSinkSource(:,:,idx) + PoincareSaddle(:,:,idx);
end
close(hWaitBar) ;

% Fixed-point categorization with analytic method
[ux,uy] = gradient(real(uv));
[vx,vy] = gradient(imag(uv));

SaddlePoincareJacobian = zeros(size(uv));
SourcePoincareJacobian = zeros(size(uv));
SinkPoincareJacobian = zeros(size(uv));
CentrePoincareJacobian = zeros(size(uv));

for idx = 1:size(uv,3)
    [row, col] = find(FixedPointsPoincare(:,:,idx) == 1);
    
    for f = 1:length(row)
        r = row(f);
        c = col(f);
        J = [ux(r,c,idx), uy(r,c,idx); vx(r,c,idx), vy(r,c,idx)];
        delta = det(J);
        tau = trace(J);
        
        [type, sp] = SourceSinkSaddle(delta, tau);
        
        if type == 2
            CentrePoincareJacobian(r,c,idx) = 1;
            
        elseif type == 1
            if sp == 1
                SourcePoincareJacobian(r,c,idx) = 1; % for spiral sources
            else
                SourcePoincareJacobian(r,c,idx) = 2; % for node sources
            end
            
        elseif type == -1
            if sp == 1
                SinkPoincareJacobian(r,c,idx) = 1;  % for spiral sinks
            else
                SinkPoincareJacobian(r,c,idx) = 2;  % for node sinks
            end
            
        elseif type == 0
            SaddlePoincareJacobian(r,c,idx) = 1;
        end
    end
end




