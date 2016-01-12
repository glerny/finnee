function [ppMS, label] = peakPickingMS_gle1(MS)
assignin('base', 'MS', MS)

%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

% Finding non nul local maxima
lmax = [1==2;  (MS(2:end-1,2) >= MS(1:end-2,2) ...
    & MS(2:end-1,2) >= MS(3:end,2)); 1 ==2] & MS(:,2) ~= 0;
% Finding local minima
lmin = [1==2;  (MS(2:end-1,2) <= MS(1:end-2,2) ...
    & MS(2:end-1,2) <= MS(3:end,2)); 1 ==2];
label = {'ind', 'MZmax', 'Imax', '-D(t)', '+D(t)',...
    'M0', 'M1', 'M2', 'Single'};

ppMS = [];
% Do peak picking
ii = 1;
while 1
    if ii > length(lmax), break, end
    if lmax(ii) && ~lmax(ii+1)
        ppMS(end+1, 1) = ii;
        
        if lmax(ii-1)
            ppMS(end, 2) = (MS(ii, 1)+ MS(ii-1,1))/2;
            ppMS(end, 3) = MS(ii, 2);
            ppMS(end, 4) =  ppMS(end, 2) - MS(ii-1, 1);
            ppMS(end, 5) = MS(ii, 1) -  ppMS(end, 2);
            
            % Find closest local minimals
            indmin1 = find(lmin(1:ii-2), 1, 'last');
            if isempty(indmin1), indmin1 = 1; end
            indmin2 = find(lmin(ii+1:end), 1, 'first')+ii;
            if isempty(indmin2), indmin2 = length(MS(:,1)); end
        else
            ppMS(end, 2) = MS(ii, 1);
            ppMS(end, 3) = MS(ii, 2);
            ppMS(end, 4) = (MS(ii, 1)-MS(ii-1, 1))/2;
            ppMS(end, 5) = (MS(ii+1, 1)-MS(ii, 1))/2;
            
            % Find closest local minimals
            indmin1 = find(lmin(1:ii-1), 1, 'last');
            if isempty(indmin1), indmin1 = 1; end
            indmin2 = find(lmin(ii+1:end), 1, 'first')+ii;
            if isempty(indmin2), indmin2 = length(MS(:,1)); end
        end
        
        ppMS(end, 6) = trapz(MS(indmin1:indmin2,1), MS(indmin1:indmin2,2));
        ppMS(end, 7) = trapz(MS(indmin1:indmin2,1), ...
            MS(indmin1:indmin2,1).*MS(indmin1:indmin2,2))/ ppMS(end, 6) ;
        ppMS(end, 8) = trapz(MS(indmin1:indmin2,1),  ...
            (MS(indmin1:indmin2,1)- ppMS(end, 7)).^2.*...
            MS(indmin1:indmin2,2))/ ppMS(end, 6) ;
        if MS(indmin1,2) == 0 &&  MS(indmin2,2) == 0
            ppMS(end, 9) = 1;
        end
        ppMS(end, 10) = indmin2 - indmin1;
        
    end
    ii = ii + 1;
end
end

