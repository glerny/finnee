function [bool, K, field] = fieldfind( structIn, label)
bool = 0; K = nan; field = {};
for ii = 1:length(structIn)
    if strcmp(structIn{ii}.label, label)
        bool = 1;
        break
    end
end
if bool
    K = ii;
    field = structIn{ii}.field;
end
end

