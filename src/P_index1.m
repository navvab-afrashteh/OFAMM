function Poincare1=P_index1(d)

D = angle(d);
s=zeros(1,4);
tap=zeros(1,4);

s(1)=D(2,1);
s(2)=D(2,2);
s(3)=D(1,2);
s(4)=D(1,1);

for i=1:4
    if i==3
        tap(i)=s(4)-s(3);
    else
        tap(i)=s(mod(i+1,4))-s(i);
    end
    if abs(tap(i))<pi/2
        tap(i)=tap(i);
    elseif tap(i)<=-pi/2
        tap(i)=tap(i)+pi;
    else
        tap(i)=tap(i)-pi;
    end
end

Poincare1=sum(tap)/pi;
