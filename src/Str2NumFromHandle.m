function Num = Str2NumFromHandle(h,DefaultNum)

str = get(h,'string');
str = str(str~=' ');
[Num, status] = str2num(str);
if ~status || Num <= 0
    Num = DefaultNum;
end

if ~strcmp(get(h,'Tag'),'TargetFramesTS')
    str = get(h,'string');
    [Num, status] = str2num(str);
    if ~status || Num <= 0
        Num = DefaultNum;
    end
end

set(h,'string',num2str(Num));