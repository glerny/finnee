function [par4cor, noise, provProf] = getBckgCorr( profileOut)
par4cor(1) = 0;
provProf = profileOut;

% if length(nonzeros(provProf(:,2))) >= 0.5*length(provProf(:,2))
%     
%     % first derivative
%     workProf = provProf;
%     workProf(:,2) = [0; (workProf(1:end-1, 2)-workProf(2:end, 2))./...
%         (workProf(1:end-1, 1)-workProf(2:end, 1))];
%     ind2null = find(workProf(:,2) == 0);
%     workProf(ind2null, :) = [];
%     indOut = find(workProf(:,2) >= mean(workProf(:,2)) + 3*std(workProf(:,2)) |...
%         workProf(:,2) <= mean(workProf(:,2)) - 3*std(workProf(:,2)));
%     while ~isempty(indOut)
%         workProf(indOut, :) = [];
%         indOut = find(workProf(:,2) >= mean(workProf(:,2)) + 3*std(workProf(:,2)) |...
%             workProf(:,2) <= mean(workProf(:,2)) - 3*std(workProf(:,2)));
%     end
%     % test null hypothesis
%     t = mean(workProf(:,2))/(std(workProf(:,2))/sqrt(length(workProf(:,2))));
%     noise(1) = std(workProf(:,2));
%     if abs(t) > 1.96
%         par4cor(1) = mean(workProf(:,2));
%         provProf(:,2) = provProf(:,2)- par4cor(1)*provProf(:,1);
%     end
%     
%     %no derivative
%     workProf = provProf;
%     ind2null = find(workProf(:,2) == 0);
%     workProf(ind2null, :) = [];
%     indOut = find(workProf >= mean(workProf) + 3*std(workProf) |...
%         workProf <= mean(workProf) - 3*std(workProf));
%     while ~isempty(indOut)
%         workProf(indOut, :) = [];
%         indOut = find(workProf >= mean(workProf) + 3*std(workProf) |...
%             workProf <= mean(workProf) - 3*std(workProf));
%     end
%     % test null hypothesis
%     t = mean(workProf)/(std(workProf)/sqrt(length(workProf)));
%     noise(3) = std(workProf);
%     if t > 1.96
%         par4cor(3) = mean(workProf);
%         provProf(:,2) = provProf(:,2)- par4cor(3);
%     end
% else
    %no derivative
    workProf = provProf;
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
    t = mean(workProf(:,2))/(std(workProf(:,2))/sqrt(length(workProf(:,2))));
    noise(1) = std(workProf(:,2));
    if t > 1.96
        par4cor(1) = mean(workProf(:,2));
        provProf(:,2) = provProf(:,2)- par4cor(1);
    end
% end
end

