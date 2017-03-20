function [SourcePoincareJacobian, SinkPoincareJacobian,...
    SourceVerifiedContours_Stack, SinkVerifiedContours_Stack] = ...
    Confirm_Contours_FPs(SourcePoincareJacobian, SinkPoincareJacobian, ...
    Contours_Stack, SourceVerifiedContours_Stack, SinkVerifiedContours_Stack,...
    Nnested_source, Nnested_sink)

[~, ~, d3] = size(SourcePoincareJacobian);
for idx = 1:d3
    C = Contours_Stack{idx};
    sourceVerified = SourceVerifiedContours_Stack{idx};
    sinkVerified = SinkVerifiedContours_Stack{idx};
    
    % source
    FPc = SourcePoincareJacobian(:,:,idx);
    [rS, cS] = find(FPc ~= 0);
    
    sourceVerified_new = [];
    for s = length(rS):-1:1
        rs = rS(s);
        cs = cS(s);
        [CorrContour] = contourConfirmation(C, sourceVerified, rs, cs, Nnested_source);
        
        if isempty(CorrContour)
            SourcePoincareJacobian(rs,cs,idx) = 0;
            continue;
        end
        sourceVerified_new = [sourceVerified_new, CorrContour];
    end
    SourceVerifiedContours_Stack{idx} = fliplr(sourceVerified_new);
    
    % sink
    FPc = SinkPoincareJacobian(:,:,idx);
    [rS, cS] = find(FPc ~= 0);
    
    sinkVerified_new = [];
    for s = length(rS):-1:1
        rs = rS(s);
        cs = cS(s);
        [CorrContour] = contourConfirmation(C, sinkVerified, rs, cs, Nnested_sink);
        if isempty(CorrContour)
            SinkPoincareJacobian(rs,cs,idx) = 0;
            continue;
        end
        sinkVerified_new = [sinkVerified_new, CorrContour];
    end
    SinkVerifiedContours_Stack{idx} = fliplr(sinkVerified_new);
end

