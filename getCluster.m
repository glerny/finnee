function results = getCluster( clustersIn, cursor_info )
%% DESCRIPTION
% 1. INTRODUCTION
% GETCLUSTER allows to obtain the different profile as well as figures of
% merits related to a single cluster. 
%
% 2. PARAMETERS:
%   .required. DOCLUSTERPLOT requires 2 parameters
%       clustersIn
%           is a structure created by the function DOCLUSTERPLOT
%       cursor_info
%           is the position in the clusters plot of the target cluster.
%           those are obtained via the 'Export Cursor Data to Workspave'
%           menu in the Data Cursor functionnalites
%
% 3. 
%
% 4. COPYRIGHT
%   Copyright 2014-2015 G. Erny (guillaume@fe.up.pt), FEUP, Porto, Portugal

%% CORE OF THE FUNCTION
% 1. INITIALISATION

fomClusters = sortrows(clustersIn.FOM.data,2);
targetInd = find(fomClusters(:,2) == cursor_info.Position(1) & ...
    fomClusters(:,6) == cursor_info.Position(2));
axeX = clustersIn.axes.time.values;
[clustP, preClustP, postClusP] = deal(zeros(length(axeX), 1));
results = [];
curCluster = clustersIn.cluster{fomClusters(targetInd, 1)}.ionicProfile;
for ii = 1:length(curCluster)
    clustP(curCluster{ii}(:,1)) = clustP(curCluster{ii}(:,1)) + ...
        curCluster{ii}(:,2)/max(curCluster{ii}(:,2));
end

if targetInd > 1
    curCluster = clustersIn.cluster{fomClusters(targetInd-1, 1)}.ionicProfile;
    for ii = 1:length(curCluster)
        preClustP(curCluster{ii}(:,1)) = preClustP(curCluster{ii}(:,1)) + ...
            curCluster{ii}(:,2)/max(curCluster{ii}(:,2));
    end
end

if targetInd < length(fomClusters(:,2))
    curCluster = clustersIn.cluster{fomClusters(targetInd+1, 1)}.ionicProfile;
    for ii = 1:length(curCluster)
        postClusP(curCluster{ii}(:,1)) = postClusP(curCluster{ii}(:,1)) + ...
            curCluster{ii}(:,2)/max(curCluster{ii}(:,2));
    end
end

figure('units','normalized','outerposition',[0 0 1 1])
iStart = find(clustP+preClustP+postClusP > 0, 1, 'first');
iStart = max(1, iStart-2);
iEnd = find(clustP+preClustP+postClusP ~= 0, 1, 'last');
iEnd = min(length(axeX), iEnd+2);

subplot(2,2,1)
plot(axeX(iStart:iEnd), clustP(iStart:iEnd), 'k')
hold on
plot(axeX(iStart:iEnd), preClustP(iStart:iEnd), 'r')
plot(axeX(iStart:iEnd), postClusP(iStart:iEnd), 'r')
hold off
title('Sum of normalized profiles in cluster')
legend('target cluster', 'pre cluster', 'post cluster')
xlabel([clustersIn.axes.time.label, ' / ',clustersIn.axes.time.unit]);
ylabel(['Normalised intensity', ' / ','profile']);

subplot(2,2,2)
MS.label = {'m/z', 'std', 'sumI', 'corr2clustP'};
MS.data = [];
curCluster = clustersIn.cluster{fomClusters(targetInd, 1)}.ionicProfile;
for ii = 1:length(curCluster)
    hold on
    prof2plot= zeros(length(axeX), 1);
    prof2plot(curCluster{ii}(:,1)) = prof2plot(curCluster{ii}(:,1)) + ...
        curCluster{ii}(:,2);
    MS.data(ii, 1) = mean(curCluster{ii}(:,3));
    MS.data(ii, 2) = std(curCluster{ii}(:,3));
    MS.data(ii, 3) = sum(curCluster{ii}(:,2));
    R = corrcoef([clustP, prof2plot]);
    MS.data(ii, 4) = R(2,1);
    plot(axeX(iStart:iEnd), prof2plot(iStart:iEnd), 'k')
end
MS.data(:, 3) = MS.data(:, 3)/max(MS.data(:, 3))*100;
hold off

title('Profiles in cluster')
xlabel([clustersIn.axes.time.label, ' / ',clustersIn.axes.time.unit]);
ylabel([clustersIn.axes.intensity.label, ' / ',clustersIn.axes.intensity.unit]);

subplot(2,2,[3 4])
stem(MS.data(:,1), MS.data(:,3), 'Marker', 'none');
title('Mass spectra')
xlabel([clustersIn.axes.mz.label, ' / ',clustersIn.axes.mz.unit]);
ylabel(['Relative intensity', ' / ',clustersIn.axes.intensity.unit]);


fs = clustersIn.prec4mz;

results.FOM.label = clustersIn.FOM.label;
results.FOM.data = fomClusters(targetInd, :);
results.MS = MS;
fprintf('\n\tTime at peak maxima: \t\t%.2f \n', fomClusters(targetInd, 2))
fprintf('\tMaximum intensity: \t%.0f \n', fomClusters(targetInd, 6))
fprintf('\tSum of max intensities: \t%.0f \n', fomClusters(targetInd, 10))
fprintf('\tSum of PIP areas: \t\t\t%.0f \n\n', fomClusters(targetInd, 11))
fprintf('\t AcuMass \t std \t Intensity \t corr2clustP\n\n')
for ii = 1:length(MS.data(:,1))
    fprintf(['\t',fs,'\t', fs, '\t%0.1f\t(%0.3f)\n'], MS.data(ii,:))
end

