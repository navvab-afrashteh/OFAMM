function [inside, Cin] = isInsideContour(c1 , c2)

x1 = c1(1,:); y1 = c1(2,:);
x2 = c2(1,:); y2 = c2(2,:);

mx = min([x1, x2]); mx = floor(mx);
Mx = max([x1, x2]); Mx = ceil(Mx);

my = min([y1, y2]); my = floor(my);
My = max([y1, y2]); My = ceil(My);

[X,Y] = meshgrid(mx:Mx, my:My);
IN1 = inpolygon(X,Y, x1,y1);
IN2 = inpolygon(X,Y, x2,y2);

area1 = sum(IN1(:));
area2 = sum(IN2(:));

minArea = min(area1, area2);

INcommon = IN1 .* IN2;

commonArea = sum(INcommon(:));

if commonArea/minArea >= 0.8
    inside = 1;
    if area1 > area2
        Cin = c2;
    else
        Cin = c1;
    end
else
    inside = 0;
    Cin = [];
end