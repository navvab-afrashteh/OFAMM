function handles = plotInstSpDir(handles, plotflag)

handles = NbinSpDir_fun(handles);
Nbins = handles.HistRose.NbinSpDir;

handles.HistRose.minSp = inf;
handles.HistRose.maxSp = 0;
if isfield(handles,'InstSpDirCLG')
    handles.InstSpDirCLG(isnan(handles.InstSpDirCLG))=0;
    handles.InstSpDirCLG(isinf(handles.InstSpDirCLG))=0;
    
    handles.HistRose.SpCLG = abs(handles.InstSpDirCLG);
    handles.HistRose.minSpCLG = quantile(handles.HistRose.SpCLG(:),0.001);
    handles.HistRose.maxSpCLG = quantile(handles.HistRose.SpCLG(:),0.999);
    handles.HistRose.AngCLG = angle(conj(handles.InstSpDirCLG));
    
    handles.HistRose.minSp = min(handles.HistRose.minSp,handles.HistRose.minSpCLG);
    handles.HistRose.maxSp = max(handles.HistRose.maxSp,handles.HistRose.maxSpCLG);
end
if isfield(handles,'InstSpDirHS')
    handles.InstSpDirHS(isnan(handles.InstSpDirHS))=0;
    handles.InstSpDirHS(isinf(handles.InstSpDirHS))=0;
    
    handles.HistRose.SpHS = abs(handles.InstSpDirHS);
    handles.HistRose.minSpHS = quantile(handles.HistRose.SpHS(:),0.001);
    handles.HistRose.maxSpHS = quantile(handles.HistRose.SpHS(:),0.999);
    handles.HistRose.AngHS = angle(conj(handles.InstSpDirHS));
    
    handles.HistRose.minSp = min(handles.HistRose.minSp,handles.HistRose.minSpHS);
    handles.HistRose.maxSp = max(handles.HistRose.maxSp,handles.HistRose.maxSpHS);
end

% calculate bin centers
binSize = (handles.HistRose.maxSp - handles.HistRose.minSp)/Nbins;
bins = (binSize/2+handles.HistRose.minSp):binSize:handles.HistRose.maxSp;

binsAll = [];
if isfield(handles,'InstSpDirCLG')
    [binsCLG, ~] = hist(handles.HistRose.SpCLG(:),bins);
    binsCLG = 100*binsCLG/length(handles.HistRose.SpCLG(:));
    binsAll = [binsAll;binsCLG];
    handles.InstSpCLG_ele = binsCLG;
    handles.InstSpCLG_bins = bins;
    
    [tCLG, rCLG] = rose(handles.HistRose.AngCLG(:),Nbins);
    rCLG = rCLG/max(rCLG);
    handles.InstDirCLG_rho = rCLG;
    handles.InstDirCLG_theta = tCLG;
end
if isfield(handles,'InstSpDirHS')
    [binsHS, ~] = hist(handles.HistRose.SpHS(:),bins);
    binsHS = 100*binsHS/length(handles.HistRose.SpHS(:));
    binsAll = [binsAll;binsHS];
    handles.InstSpHS_ele = binsHS;
    handles.InstSpHS_bins = bins;
    
    [tHS, rHS] = rose(handles.HistRose.AngHS(:),Nbins);
    rHS = rHS/max(rHS);
    handles.InstDirHS_rho = rHS;
    handles.InstDirHS_theta = tHS;
end

% plotting
if plotflag && (isfield(handles,'InstSpDirCLG') || isfield(handles,'InstSpDirHS'))
    if isfield(handles,'InstSpDirFig')
        if ishandle(handles.InstSpDirFig)
            figure(handles.InstSpDirFig);
        else
            handles.InstSpDirFig = figure(102);
            set(handles.InstSpDirFig,'visible','off');
            set(handles.InstSpDirFig,'numbertitle','off','name','Instantaneous Speeds and Angles histogram for selected ROI','units','inches','position',[1 3.5 6 2.5]);
            set(handles.InstSpDirFig,'visible','on');
        end
    else
        handles.InstSpDirFig = figure(102);
        set(handles.InstSpDirFig,'visible','off');
        set(handles.InstSpDirFig,'numbertitle','off','name','Instantaneous Speeds and Angles histogram for selected ROI','units','inches','position',[1 3.5 6 2.5]);
        set(handles.InstSpDirFig,'visible','on');
    end
    c1 = [0.08,0.17,0.55];
    c2 = [0.55,0,0];
    figure(handles.InstSpDirFig);
    subplot(121); cla;
    h_bar = bar(bins,binsAll');
    if length(h_bar)>1
        set(h_bar(1),'facecolor',c1)
        set(h_bar(2),'facecolor',c2)
    elseif length(h_bar)==1
        if isfield(handles,'InstSpDirCLG')
            set(h_bar,'facecolor',c1)
        else
            set(h_bar,'facecolor',c2)
        end
    end
    
    maxY = max(binsAll(:));
    xlim([min(bins)-binSize, max(bins)+binSize]);
    ylim([0 maxY*1.1]);
    set(gca,'box','off','TickDir','out');
    xlabel(sprintf('Inst. Speed (p/f)')); ylabel('Percentage');
    
    figure(handles.InstSpDirFig);
    subplot(122); cla;
    if isfield(handles,'InstSpDirCLG')
        h = polar(tCLG,rCLG);
        set(h,'color',c1);
        hold on;
        if isfield(handles,'InstSpDirHS')
            h = polar(tHS,rHS);
            set(h,'color',c2);
            hold off;
            subplot(121);
            legend('CLG','HS','Location','best')
        else
            subplot(121);
            legend('CLG','Location','best')
        end
    else
        if isfield(handles,'InstSpDirHS')
            h = polar(tHS,rHS);
            set(h,'color',c2);
            hold off;
            subplot(121);
            legend('HS','Location','best')
        end
    end
end
