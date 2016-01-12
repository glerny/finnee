function [ noise, signal ] = getNoise_s(data)
nbrPtMin = 4;
Dstd = 0.05;
signal = max(data);
conv = false;
Ddl = abs(data(1:end-1)- data(2:end));
std1 = std(Ddl);

while 1
    if length(Ddl) <= nbrPtMin
        conv = false;
        break
    end
    ind = find(Ddl > (max(Ddl) + mean(Ddl))/2);
    Ddl(ind) = [];
    
    if abs((std1-std(Ddl))/std1) < Dstd
        conv = true;
        break
    else
        std1 = std(Ddl);
    end
end

if conv
    noise = max(Ddl);
else
    noise = nan;
end
    
