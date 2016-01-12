function [par4cor, noise, profileOut] = getBckgCorrNew( profileIn, m, n)
% m % of non zeros values to use polynomial otherwise cst
% n % degree of polynomial fit

profileOut = profileIn;
workProf = profileOut;

if length(nonzeros(profileOut(:,2))) >= m*length(profileOut(:,2)) %use polyfit with degree n
    ind2null = find(workProf(:,2) == 0);
    workProf(ind2null, :) = [];
    p = polyfit(workProf(:,1), workProf(:,2), n);
    profileRsd = workProf(:,2) - polyval(p, workProf(:,1));
    indOut = find(profileRsd >= mean(profileRsd) + 3*std(profileRsd) |...
        profileRsd <= mean(profileRsd) - 3*std(profileRsd));
    while ~isempty(indOut) && length(profileRsd) >= m/2*length(profileOut(:,2)) 
        workProf(indOut, :) = [];
        p = polyfit(workProf(:,1), workProf(:,2), n);
    	profileRsd = workProf(:,2) - polyval(p, workProf(:,1));
    	indOut = find(profileRsd >= mean(profileRsd) + 3*std(profileRsd) |...
     	   profileRsd <= mean(profileRsd) - 3*std(profileRsd));
    end
    noise = std(profileRsd);
    par4cor = p;
    profileOut(:,2) = profileOut(:,2)- polyval(p, profileOut(:,1));
    
else
    %use constant
    ind2null = find(workProf(:,2) == 0);
    workProf(ind2null, :) = [];
    indOut = find(workProf(:,2) >= mean(workProf(:,2)) + 3*std(workProf(:,2)) |...
        workProf(:,2) <= mean(workProf(:,2)) - 3*std(workProf(:,2)));
    while ~isempty(indOut)
        workProf(indOut, :) = [];
        indOut = find(workProf(:,2) >= mean(workProf(:,2)) + 3*std(workProf(:,2)) |...
	   workProf(:,2) <= mean(workProf(:,2)) - 3*std(workProf(:,2)));
    end

    % test null hypothesis
    t = mean(workProf)/(std(workProf)/sqrt(length(workProf)));
    noise = std(workProf(:,2));
    par4cor(n+1) = mean(workProf(:,2));
    profileOut(:,2) = profileOut(:,2)- par4cor(n+1);
end
end

