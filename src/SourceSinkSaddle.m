function [type, sp] = SourceSinkSaddle(delta, tau)

% This function is to identify if a fixed-point is source, sink, or saddle
% point. 'delta' is the determinant of Jcobian matrix,J, and 'tau' is the 
% trace of J. If delta is less than zero, the point is saddle point. If 
% tau > 0 and delta > 0, then the point is a source. If tau*tau < 4*delta, 
% then the source is a spiral source. If not, this is an ordinary source. 
% On the other hand, if tau < 0 and delta > 0, the point is a sink. Again,
% if tau*tau < 4*delta, then the sink is a spiral one. If not, this is an 
% ordinary sink. 'Centres' are fixed-points which tau is zero and delta is
% greater than zero.

if delta < 0
    
    type = 0;
    sp = 0;
    
elseif delta > 0
    
    if tau > 0
        type = 1;
        if tau*tau < 4*delta
            sp = 1;
        else
            sp = 0;
        end
        
    elseif tau < 0
        type = -1;
        if tau*tau < 4*delta
            sp = 1;
        else
            sp = 0;
        end
    elseif tau == 0
        type = 2;
        sp = 1;
    end
    
elseif delta == 0
    type = NaN;
    sp = NaN;
end