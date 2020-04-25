function handles = plotTempSpTrajLen(handles, plotflag)

handles = NbinSpTraj_fun(handles);
Nbins = handles.TrajHist.NbinSpTraj;

handles.TrajHist.minSp = inf;
handles.TrajHist.maxSp = 0;
handles.TrajHist.minLen = inf;
handles.TrajHist.maxLen = 0;
handles.TrajHist.minMaxSp = inf;
handles.TrajHist.maxMaxSp = 0;
if isfield(handles,'TempSpTrajLenCLG')
    AllSpCLG = zeros(length(handles.TempSpTrajLenCLG.speed{1}),length(handles.TempSpTrajLenCLG.speed));
    totalDistCLG = zeros(1,length(handles.TempSpTrajLenCLG.speed));
    for point = 1:length(handles.TempSpTrajLenCLG.speed)
        AllSpCLG(:,point) = handles.TempSpTrajLenCLG.speed{point};
        
        str = handles.TempSpTrajLenCLG.str{point};
        xs = str(:,1); ys = str(:,2);
        vects = xs + 1i*ys;
        diffVects = diff(vects);
        distances = abs(diffVects);
        totalDistCLG(point) = sum(distances);
    end
    
    tSpCLG = handles.TempSpTrajLenCLG.tspeed{point};
    AllSpCLG(isnan(AllSpCLG))=0;
    AllSpCLG(isinf(AllSpCLG))=0;
    totalDistCLG(isnan(totalDistCLG))=0;
    totalDistCLG(isinf(totalDistCLG))=0;
    AllMaxSpCLG = max(AllSpCLG);
    
    handles.TrajHist.minSpCLG = quantile(AllSpCLG(:),0.001);
    handles.TrajHist.maxSpCLG = quantile(AllSpCLG(:),0.999);
    handles.TrajHist.minSp = min(handles.TrajHist.minSp,handles.TrajHist.minSpCLG);
    handles.TrajHist.maxSp = max(handles.TrajHist.maxSp,handles.TrajHist.maxSpCLG);
    
    handles.TrajHist.minLenCLG = quantile(totalDistCLG(:),0.001);
    handles.TrajHist.maxLenCLG = quantile(totalDistCLG(:),0.999);
    handles.TrajHist.minLen = min(handles.TrajHist.minLen,handles.TrajHist.minLenCLG);
    handles.TrajHist.maxLen = max(handles.TrajHist.maxLen,handles.TrajHist.maxLenCLG);
    
    handles.TrajHist.minMaxSpCLG = quantile(AllMaxSpCLG,0.001);
    handles.TrajHist.maxMaxSpCLG = quantile(AllMaxSpCLG,0.999);
    handles.TrajHist.minMaxSp = min(handles.TrajHist.minMaxSp,handles.TrajHist.minMaxSpCLG);
    handles.TrajHist.maxMaxSp = max(handles.TrajHist.maxMaxSp,handles.TrajHist.maxMaxSpCLG);
end

if isfield(handles,'TempSpTrajLenHS')
    AllSpHS = zeros(length(handles.TempSpTrajLenHS.speed{1}),length(handles.TempSpTrajLenHS.speed));
    totalDistHS = zeros(1,length(handles.TempSpTrajLenHS.speed));
    for point = 1:length(handles.TempSpTrajLenHS.speed)
        AllSpHS(:,point) = handles.TempSpTrajLenHS.speed{point};
        
        str = handles.TempSpTrajLenHS.str{point};
        xs = str(:,1); ys = str(:,2);
        vects = xs + 1i*ys;
        diffVects = diff(vects);
        distances = abs(diffVects);
        totalDistHS(point) = sum(distances);
    end
    AllMaxSpHS = max(AllSpHS);
    tSpHS = handles.TempSpTrajLenHS.tspeed{point};
    AllSpHS(isnan(AllSpHS))=0;
    AllSpHS(isinf(AllSpHS))=0;
    totalDistHS(isnan(totalDistHS))=0;
    totalDistHS(isinf(totalDistHS))=0;
    
    handles.TrajHist.minSpHS = quantile(AllSpHS(:),0.001);
    handles.TrajHist.maxSpHS = quantile(AllSpHS(:),0.999);
    handles.TrajHist.minSp = min(handles.TrajHist.minSp,handles.TrajHist.minSpHS);
    handles.TrajHist.maxSp = max(handles.TrajHist.maxSp,handles.TrajHist.maxSpHS);
    
    handles.TrajHist.minLenHS = quantile(totalDistHS(:),0.001);
    handles.TrajHist.maxLenHS = quantile(totalDistHS(:),0.999);
    handles.TrajHist.minLen = min(handles.TrajHist.minLen,handles.TrajHist.minLenHS);
    handles.TrajHist.maxLen = max(handles.TrajHist.maxLen,handles.TrajHist.maxLenHS);
    
    handles.TrajHist.minMaxSpHS = quantile(AllMaxSpHS,0.001);
    handles.TrajHist.maxMaxSpHS = quantile(AllMaxSpHS,0.999);
    handles.TrajHist.minMaxSp = min(handles.TrajHist.minMaxSp,handles.TrajHist.minMaxSpHS);
    handles.TrajHist.maxMaxSp = max(handles.TrajHist.maxMaxSp,handles.TrajHist.maxMaxSpHS);
end

% calculate bin centers
binSizeSp = (handles.TrajHist.maxSp - handles.TrajHist.minSp)/Nbins;
binsSp = (binSizeSp/2+handles.TrajHist.minSp):binSizeSp:handles.TrajHist.maxSp;

binSizeLen = (handles.TrajHist.maxLen - handles.TrajHist.minLen)/Nbins;
binsLen = (binSizeLen/2+handles.TrajHist.minLen):binSizeLen:handles.TrajHist.maxLen;

binSizeMaxSp = (handles.TrajHist.maxMaxSp - handles.TrajHist.minMaxSp)/Nbins;
binsMaxSp = (binSizeMaxSp/2+handles.TrajHist.minMaxSp):binSizeMaxSp:handles.TrajHist.maxMaxSp;

binsAllSp = [];
binsAllLen = [];
binsAllMaxSp = [];
if isfield(handles,'TempSpTrajLenCLG')
    [binsSpCLG, ~] = hist(AllSpCLG(:),binsSp);
    binsSpCLG = 100*binsSpCLG/length(AllSpCLG(:));
    binsAllSp = [binsAllSp;binsSpCLG];
    handles.TempSpCLG_ele = binsSpCLG;
    handles.TempSpCLG_bins = binsSp;
    
    [binsLenCLG, ~] = hist(totalDistCLG,binsLen);
    binsLenCLG = 100*binsLenCLG/length(totalDistCLG);
    binsAllLen = [binsAllLen;binsLenCLG];
    handles.TempLenCLG_ele = binsLenCLG;
    handles.TempLenCLG_bins = binsLen;
    
    [binsMaxCLG, ~] = hist(AllMaxSpCLG(:),binsMaxSp);
    binsMaxCLG = 100*binsMaxCLG/length(AllMaxSpCLG(:));
    binsAllMaxSp = [binsAllMaxSp;binsMaxCLG];
    handles.MaxTempSpCLG_ele = binsMaxCLG;
    handles.MaxTempSpCLG_bins = binsMaxSp;
end

if isfield(handles,'TempSpTrajLenHS')
    [binsSpHS, ~] = hist(AllSpHS(:),binsSp);
    binsSpHS = 100*binsSpHS/length(AllSpHS(:));
    binsAllSp = [binsAllSp;binsSpHS];
    handles.TempSpHS_ele = binsSpHS;
    handles.TempSpHS_bins = binsSp;
    
    [binsLenHS, ~] = hist(totalDistHS,binsLen);
    binsLenHS = 100*binsLenHS/length(totalDistHS);
    binsAllLen = [binsAllLen;binsLenHS];
    handles.TempLenHS_ele = binsLenHS;
    handles.TempLenHS_bins = binsLen;
    
    [binsMaxHS, ~] = hist(AllMaxSpHS(:),binsMaxSp);
    binsMaxHS = 100*binsMaxHS/length(AllMaxSpHS(:));
    binsAllMaxSp = [binsAllMaxSp;binsMaxHS];
    handles.MaxTempSpHS_ele = binsMaxHS;
    handles.MaxTempSpHS_bins = binsMaxSp;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plotting
if plotflag && (isfield(handles,'TempSpTrajLenCLG') || isfield(handles,'TempSpTrajLenHS'))
    if isfield(handles,'TempSptrajLenFig')
        if ishandle(handles.TempSptrajLenFig)
            figure(handles.TempSptrajLenFig);
        else
            handles.TempSptrajLenFig = figure(103);
            set(handles.TempSptrajLenFig,'visible','off');
            set(handles.TempSptrajLenFig,'numbertitle','off','name','Temporal Speeds, profile and histogram, and Trajectory Length histogram for selected ROI','units','inches','position',[1 1 12 2.5]);
            set(handles.TempSptrajLenFig,'visible','on');
        end
    else
        handles.TempSptrajLenFig = figure(103);
        set(handles.TempSptrajLenFig,'visible','off');
        set(handles.TempSptrajLenFig,'numbertitle','off','name','Temporal Speeds, profile and histogram, and Trajectory Length histogram for selected ROI','units','inches','position',[1 1 12 2.5]);
        set(handles.TempSptrajLenFig,'visible','on');
    end
    clf
    % temporal profile
    c1 = [0.08,0.17,0.55];
    c2 = [0.55,0,0];
    figure(handles.TempSptrajLenFig);
    subplot(141); cla;
    if isfield(handles,'TempSpTrajLenCLG')
        LegendCLG = plot(tSpCLG,AllSpCLG,'color',c1);hold on;
        if isfield(handles,'TempSpTrajLenHS')
            LegendHS = plot(tSpHS,AllSpHS,'color',c2);
            legend([LegendCLG(1),LegendHS(1)],'CLG','HS','Location','best')
            hold off;
        else
            legend('CLG','Location','best')
        end
    else
        if isfield(handles,'TempSpTrajLenHS')
            plot(tSpHS,AllSpHS,'color',c2);
            legend('HS','Location','best')
            hold off;
        end
    end
    ylabel('Temp. Speed (p/f)');
    xlabel('Frame Number');
    
    % teporal speed histogram
    figure(handles.TempSptrajLenFig);
    subplot(142); cla;
    
    h_bar = bar(binsSp,binsAllSp');
    if length(h_bar)>1
        set(h_bar(1),'facecolor',c1)
        set(h_bar(2),'facecolor',c2)
    elseif length(h_bar)==1
        if isfield(handles,'TempSpTrajLenCLG')
            set(h_bar,'facecolor',c1)
        else
            set(h_bar,'facecolor',c2)
        end
    end
    
    maxY = max(binsAllSp(:));
    xlim([min(binsSp)-binSizeSp, max(binsSp)+binSizeSp]);
    ylim([0 maxY*1.1]);
    set(gca,'box','off','TickDir','out');
    xlabel(sprintf('Temp. Speed (p/f)')); ylabel('Percentage');
    
    % Maximum teporal speed histogram
    figure(handles.TempSptrajLenFig);
    subplot(143); cla;
    
    h_bar = bar(binsMaxSp,binsAllMaxSp');
    if length(h_bar)>1
        set(h_bar(1),'facecolor',c1)
        set(h_bar(2),'facecolor',c2)
    elseif length(h_bar)==1
        if isfield(handles,'TempSpTrajLenCLG')
            set(h_bar,'facecolor',c1)
        else
            set(h_bar,'facecolor',c2)
        end
    end
    
    maxY = max(binsAllMaxSp(:));
    xlim([min(binsMaxSp)-binSizeMaxSp, max(binsMaxSp)+binSizeMaxSp]);
    ylim([0 maxY*1.1]);
    set(gca,'box','off','TickDir','out');
    xlabel(sprintf('Max Temp. Speed (p/f)')); ylabel('Percentage');
    
    % Trajectory Distance histogram
    figure(handles.TempSptrajLenFig);
    subplot(144); cla;
    
    h_bar = bar(binsLen,binsAllLen');
    if length(h_bar)>1
        set(h_bar(1),'facecolor',c1)
        set(h_bar(2),'facecolor',c2)
    elseif length(h_bar)==1
        if isfield(handles,'TempSpTrajLenCLG')
            set(h_bar,'facecolor',c1)
        else
            set(h_bar,'facecolor',c2)
        end
    end
    
    maxY = max(binsAllLen(:));
    xlim([min(binsLen)-binSizeLen, max(binsLen)+binSizeLen]);
    ylim([0 maxY*1.1]);
    set(gca,'box','off','TickDir','out');
    xlabel(sprintf('Trajectory Length (pixels)')); ylabel('Percentage');
end
