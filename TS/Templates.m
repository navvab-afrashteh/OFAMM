function [T1, T2, T3, T4, dir] = Templates(k, points)

len = length(points);
direction = zeros(len,1);

for idx= 1:len
    direction(idx) = points(idx,2) + 1i*points(idx,1);
end
direction = direction - (k+1)*(1+1i);

T1 = zeros(2*len,1);
T2 = zeros(2*len,1);

T1(1:len) = real(direction);
T2(1:len) = imag(direction);

ring = circshift(direction, [-1, 0]) - direction;

T1(len+1:2*len) = real(ring);
T2(len+1:2*len) = imag(ring);

T3 = [abs(direction); zeros(len,1)];
T4 = [zeros(len,1); abs(ring)];

dir = [direction ; ring ];

T1 = 1 ./ T1;
T2 = 1 ./ T2;
T3 = 1 ./ T3;
T4 = 1 ./ T4;

T1inf = isinf(T1);
T2inf = isinf(T2);
T3inf = isinf(T3);
T4inf = isinf(T4);

for idx = 1:2*len
    if T1inf(idx) == 1
        T1(idx) = 0;
    end
    if T2inf(idx) == 1
        T2(idx) = 0;
    end
    if T3inf(idx) == 1
        T3(idx) = 0;
    end
    if T4inf(idx) == 1
        T4(idx) = 0;
    end
end
