function [source, sink, C] = refiningContours(u,v,Mask,Nmin, Lmin_source, Lmax_sink)

[ux, ~] = gradient(u);
[~, vy] = gradient(v);

div1 = ux + vy; % Divergence
div1 = div1 .* Mask;

[C,h] = contour(div1);
set(h,'visible','off')

lCS = size(C, 2);

iL = [];
L = [];
k = 1;
while (k < lCS)
    iL = [iL, k];
    k = k + C(2, k) + 1;
end

Nsource=0;
source = [];
for idx = 1:length(iL)
    k = iL(idx);
    level = C(1,k);
    N = C(2,k);
    
    x = C(1,k+1:k+N);
    y = C(2,k+1:k+N);
    if N > Nmin && level > Lmin_source && x(1) == x(end) && y(1) == y(end)%&& Ic > Imean %&& level < 1.5
        Nsource=Nsource+1;
        source = [source, [k; N; level]];
    end
end

Nsink=0;
sink =[];
for idx = 1:length(iL)
    k = iL(idx);
    level = C(1,k);
    N = C(2,k);
    
    x = C(1,k+1:k+N);
    y = C(2,k+1:k+N);
    if N > Nmin && level < Lmax_sink && x(1) == x(end) && y(1) == y(end) %&& Ic > Imean %&& level > -1.5
        Nsink = Nsink+1;
        sink = [sink, [k; N; level]];
    end
end

