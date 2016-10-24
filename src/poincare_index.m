function [SinkSource,Saddle]=poincare_index(uv)

[row, col]=size(uv);

SinkSource=zeros(size(uv));
Saddle=zeros(size(uv));

PoincareIdx = nlfilter(uv,[2 2],@P_index1);

for i=1:row
    for j=1:col
        if PoincareIdx(i,j)> 0.9 
            SinkSource(i,j)=1;
        elseif PoincareIdx(i,j)< -0.9
            Saddle(i,j)=1;
        end
    end
end

