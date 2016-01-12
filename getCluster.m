function results = getCluster( dataOut, cursor_info )
%% DESCRIPTION
% 1. INTRODUCTION
% GETCLUSTER allows to obtain the different profile as well as figures of
% merits related to a single cluster. 
%
% 2. PARAMETERS:
%   .required. DOCLUSTERPLOT requires 2 parameters
%       dataOut
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
fidReadDat = fopen(dataOut.results.path2DatFile, 'rb');
fomClusters.label = dataOut.results.fomClusters.label;
fomClusters.data = sortrows(dataOut.results.fomClusters.data,2);
targetInd = find(fomClusters.data(:,2) == cursor_info.Position(1) & ...
    fomClusters.data(:,6) == cursor_info.Position(2));
axeX = dataOut.results.axeX;

figure('units','normalized','outerposition',[0 0 1 1])
subplot(2,2,1)
if targetInd-1 >= 1
    hold on
    listPIP = dataOut.results.selClusters{fomClusters.data(targetInd-1,1)};
    data2plot = zeros(length(axeX), 1);
    index = dataOut.results.index2DotDat(listPIP(1), :);
    fseek(fidReadDat, index(1), 'bof');
    PIP = fread(fidReadDat, [(index(2)-index(1))/(index(3)*8), ...
        index(3)], 'double');
    data2plot(PIP(:,1)) = PIP(:,2);
    [~, iM] = max (data2plot);
    iS = 1; %max(1, find(data2plot(1:iM)==0, 1, 'last') - 1);
    iE = length(data2plot); %min(length(data2plot), find(data2plot(iM:end)==0, 1, 'first') + iM);
    plot(axeX(iS:iE),data2plot(iS:iE), 'r');
end
if targetInd+1 <= length(fomClusters.data(:,1))
    hold on
    listPIP = dataOut.results.selClusters{fomClusters.data(targetInd+1,1)};
    data2plot = zeros(length(axeX), 1);
    index = dataOut.results.index2DotDat(listPIP(1), :);
    fseek(fidReadDat, index(1), 'bof');
    PIP = fread(fidReadDat, [(index(2)-index(1))/(index(3)*8), ...
        index(3)], 'double');
    data2plot(PIP(:,1)) = PIP(:,2);
    [~, iM] = max (data2plot);
    iS = 1; %max(1, find(data2plot(1:iM)==0, 1, 'last') - 1);
    iE = length(data2plot); %min(length(data2plot), find(data2plot(iM:end)==0, 1, 'first') + iM);
    plot(axeX(iS:iE),data2plot(iS:iE), 'r');
end
hold on
listPIP = dataOut.results.selClusters{fomClusters.data(targetInd,1)};
data2plot = zeros(length(axeX), 1);
index = dataOut.results.index2DotDat(listPIP(1), :);
fseek(fidReadDat, index(1), 'bof');
PIP = fread(fidReadDat, [(index(2)-index(1))/(index(3)*8), ...
    index(3)], 'double');
data2plot(PIP(:,1)) = PIP(:,2);
[~, iM] = max (data2plot);
iS = 1; %max(1, find(data2plot(1:iM)==0, 1, 'last') - 1);
iE = length(data2plot); %min(length(data2plot), find(data2plot(iM:end)==0, 1, 'first') + iM);
plot(axeX(iS:iE),data2plot(iS:iE), 'k');

subplot(2,2,2)
MSSpectra =[];
text2disp = {};
fS = dataOut.results.formatSpec;
listofPip = [];
for ii = 1:length(listPIP)
   hold on
   data2plot = zeros(length(axeX), 1);
   index = dataOut.results.index2DotDat(listPIP(ii), :);
   fseek(fidReadDat, index(1), 'bof');
   PIP = fread(fidReadDat, [(index(2)-index(1))/(index(3)*8), ...
       index(3)], 'double');
   data2plot(PIP(:,1)) = PIP(:,2);
   [~, iM] = max (data2plot);
      iS = 1; %max(1, find(data2plot(1:iM)==0, 1, 'last') - 1);
    iE = length(data2plot); %min(length(data2plot), find(data2plot(iM:end)==0, 1, 'first') + iM);
    plot(axeX(iS:iE),data2plot(iS:iE), 'k');
   MSSpectra(end+1, 1) = dataOut.results.fom.data(listPIP(ii), 2);
   MSSpectra(end, 2) = dataOut.results.fom.data(listPIP(ii), 5);
   text2disp{ii} = [num2str(dataOut.results.fom.data(listPIP(ii), 2), fS), ' \t '... 
       '(', num2str(dataOut.results.fom.data(listPIP(ii), 3), fS),')',' \t ', ...
       num2str(dataOut.results.fom.data(listPIP(ii), 5), '%.0f')];
   assignin('base', 'data2plot', data2plot)
   assignin('base', 'listofPip', listofPip)
   listofPip(:,end+1) = data2plot/max(data2plot);
end
R = corrcoef(listofPip);
mean(mean(R))
subplot(2,2,[3 4])
stem(MSSpectra(:,1), MSSpectra(:,2), 'Marker', 'none');
assignin('base', 'listofPip', listofPip)

fprintf('\n\tTime at peak maxima: \t\t%.2f \n', fomClusters.data(targetInd, 2))
fprintf('\tSum of max intensities: \t%.0f \n', fomClusters.data(targetInd, 10))
fprintf('\tSum of PIP areas: \t\t\t%.0f \n\n', fomClusters.data(targetInd, 11))
fprintf('\t AcuMass \t (std) \t Intensity \n\n')
for ii = 1:length(text2disp)
    fprintf(['\t', text2disp{ii}, '\n'])
end

results.listPIP = listPIP;
results.text2disp = text2disp;

