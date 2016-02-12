function clustersOut = doClustersPlot(finneeStc, address, HL, minInt, varargin)
%% DESCRIPTION
% 1. INTRODUCTION
% DOCLUSTERPLOT allows to select clusters at given hierarchical level
% defined as the minimal correlation coefficient between pure ion profiles
% (PIP) within each cluster.
%
% 2. INPUT PARAMETERS
%   .required. DOCLUSTERPLOT requires at least 4 parameters
%       finneeStc
%           is the finnee structure that contain information about the run
%           and link and indexation of the associated dat file. The
%           strcuture should have been create by function such as 
%           MZML2STRUCT
%       address
%           defined the target cluster analysis (HACA) within the target 
%           dataset. The format should be 'HA@dataset'.
%       HL
%           is the hierarchical level at which the clusters will be
%           retrieved from the previous HACA analysis (see DOHACA). The
%           value should be beetween  1 (all PIP in individual clusters)
%           and 0 (all PIP in one cluster). The recomand value is beetween
%           0.95 and 0.9.
%       minInt
%           is the threshold intensity. Clusters that do not contain one
%           PIP with a maximum intensity higher than this threshold will
%           not be displayed
%
%   .optionals. VARARGIN describes the optional paramters.  
%       'displayOff' 
%           Allow to avoid any display in the matlab command window
%       'minPIPperCluster' followed by an integer (default 2)
%           Only record clusters with 'minPIPperCluster' PIP
%
% 3. OUTPUT PARAMETERS
%
% 4. EXAMPLES:
%	dataOut = doClustersPlot(finneeStc, '1@3', 0.95, 1000)
%
% 5. COPYRIGHT
%   Copyright 2015-2016 G. Erny (guillaume@fe.up.pt), FEUP, Porto, Portugal

%% CORE OF THE FUNCTION
% 1. INITIALISATION
info.function.functionName =  'doHACA';
info.function.description{1} = 'Calculate the herarchical structure of '' dataset';
info.function.matlabVersion = '8.5.0.197613 (R2015a)';
info.function.version = '18/01/2016';
info.function.ownerContact = 'guillaume@fe.up.pt';

[parameters, options] = ...
    initFunction(nargin, finneeStc, address, HL, minInt, varargin);
%INITFUNCTION - sub function used to verify the entries and load the optional and
% compulsory and optional parameters

m = parameters.dataset;
n = parameters.HStruct;

fidReadDat = fopen(finneeStc.path2dat, 'rb');
fom.label = finneeStc.dataset{m}.FOM.label;
fom.data = finneeStc.dataset{m}.FOM.data;
axeX = finneeStc.dataset{m}.axes.time.values;
        
formatSpec = finneeStc.info.parameters.prec4mz;

clusters.step = finneeStc.dataset{m}.CA{n}.HStructure.step;

% 2. FIND THE HL AND LOAD TARGETED CLUSTERS
targetStep = [];
for ii = 1:length(clusters.step)
    if clusters.step{ii}.HL < parameters.HL && isempty(targetStep)
        targetStep = ii-1;
    end
end
if isempty(targetStep), targetStep = ii ; end


index = clusters.step{targetStep}.index2DotDat;
fseek(fidReadDat, index(1), 'bof');
data2load = fread(fidReadDat, [(index(2)-index(1))/(index(3)*8), ...
    index(3)], 'double');
selClusters{length(data2load(:,1))} = {};
for ii = 1:length(data2load(:,1))
    selClusters{ii} = nonzeros(data2load(ii,:));
end

% 3. FILTERING AND CALCULATING CLUSTERS FIGURE OF MERITS
[sumPIP, sumCluSIP] = deal(zeros(length(axeX), 1));
fomClusters.label = {'clusterID' 'meanTmax' 'std' 'meanM1' 'std'...
    'massSIPmax' 'std' 'MassCentroid' 'IatSIPMax' 'sumI' 'sumArea'} ;
fomClusters.data = [];

clustersOut.title = ['Clusters @ HL = ', num2str(HL)] ;
clustersOut.axes = finneeStc.dataset{m}.axes;
clustersOut.prec4mz = finneeStc.info.parameters.prec4mz;
clustersOut.cluster = {};
for ii = 1:length(selClusters)
    clustersOut.cluster{ii}.ionicProfile = {};
    sumPIPinCl = zeros(length(axeX), 1);
    maxInt = 0;
    for jj = 1:length(selClusters{ii})
        index = finneeStc.dataset{m}.indexInDat(selClusters{ii}(jj), :);
        fseek(fidReadDat, index(1), 'bof');
        PIP = fread(fidReadDat, [(index(2)-index(1))/(index(3)*8), ...
            index(3)], 'double');
        clustersOut.cluster{ii}.ionicProfile{jj} = PIP;
        sumPIPinCl(PIP(:,1)) = sumPIPinCl(PIP(:,1)) + PIP(:,2);
        if maxInt < max(PIP(:,2)), maxInt =  max(PIP(:,2)); end
    end
    
    if maxInt >= parameters.minIntensity
        sumCluSIP = sumCluSIP + sumPIPinCl;
        fomClusters.data(end+1, 1) = ii;
        fomClusters.data(end, 2) = mean(fom.data(selClusters{ii}, 6));
        fomClusters.data(end, 3) = std(fom.data(selClusters{ii}, 6));
        fomClusters.data(end, 4) = mean(fom.data(selClusters{ii}, 8));
        fomClusters.data(end, 5) = std(fom.data(selClusters{ii}, 8));
        fomClusters.data(end, 6) = fom.data(selClusters{ii}(1), 2);
        fomClusters.data(end, 7) = fom.data(selClusters{ii}(1), 3);
        fomClusters.data(end, 8) = sum(fom.data(selClusters{ii}, 2).*...
            fom.data(selClusters{ii}, 5))/sum(fom.data(selClusters{ii}, 5));
        fomClusters.data(end, 9) = fom.data(selClusters{ii}(1), 5);
        fomClusters.data(end, 10) = sum(fom.data(selClusters{ii}, 5));
        fomClusters.data(end, 11) = sum(fom.data(selClusters{ii}, 7));
    end
    sumPIP = sumPIP + sumPIPinCl;
end
fclose(fidReadDat);
clustersOut.FOM = fomClusters;

timeLabel = finneeStc.dataset{m}.axes.time.label;
timeUnit = finneeStc.dataset{m}.axes.time.unit;
mzLabel = finneeStc.dataset{m}.axes.mz.label;
mzUnit = finneeStc.dataset{m}.axes.mz.unit;
intLabel = finneeStc.dataset{m}.axes.intensity.label;
intUnit = finneeStc.dataset{m}.axes.intensity.unit;

clustersOut.trace{1}.title = 'Sum final profiles in clusters';
clustersOut.trace{1}.data = [axeX sumPIP];
clustersOut.trace{1}.dateOfCreation = datetime;
clustersOut.trace{1}.plotType = 'profile';
clustersOut.trace{1}.axeX.label = timeLabel;
clustersOut.trace{1}.axeX.unit = timeUnit;
clustersOut.trace{1}.axeY.label = intLabel;
clustersOut.trace{1}.axeY.unit = intUnit;

clustersOut.trace{2}.title =['Sum profiles in clusters with intensity >'...
    num2str(minInt)];
clustersOut.trace{2}.data = [axeX sumCluSIP];
clustersOut.trace{2}.dateOfCreation = datetime;
clustersOut.trace{2}.plotType = 'profile';
clustersOut.trace{2}.axeX.label = timeLabel;
clustersOut.trace{2}.axeX.unit = timeUnit;
clustersOut.trace{2}.axeY.label = intLabel;
clustersOut.trace{2}.axeY.unit = intUnit;

% 4. CLUSTERS PLOT
listMZ2plot3D = fomClusters.data(:,6);
listI2plot3D =  fomClusters.data(:,9);
listT2plot3D = fomClusters.data(:,2);
figure('units','normalized','outerposition',[0 0 1 1], 'Name', ... % this does make it full size on my computer don't knwo why
    ['Clusters plots of ',finneeStc.info.parameters.fileID, ' @ ',...
    address, ' with HL = ', num2str(parameters.HL), ...
    ' and intenisty of clusters >= ', num2str(parameters.minIntensity)])

subplot(5,1,1)
plot(axeX, sumPIP, 'k')
hold on
plot(axeX, sumCluSIP, 'r')
hold off
legend('all clusters', ['clust. int. > ', num2str(minInt)])  
title('Sum of profiles in clusters');
xlabel([timeLabel,' / ', timeUnit]);
ylabel([intLabel,' / ', intUnit]);
v1 = axis;

subplot(5,1, 2:5) % CLUSTERS PLOT
codeSize = int32(listI2plot3D/max(listI2plot3D)*(500) + 5);
hold on
scatter(listT2plot3D, listMZ2plot3D, ...
    codeSize, 'k');
title('Clusters Plot');
xlabel([timeLabel,' / ', timeUnit]);
ylabel([mzLabel,' / ', mzUnit]);
v2 = axis;
axis([v1(1) v1(2) v2(3) v2(4)])
hold off

clustersOut.plot.x.value = listT2plot3D;
clustersOut.plot.x.label = [timeLabel,' / ', timeUnit];
clustersOut.plot.y.value = listMZ2plot3D;
clustersOut.plot.y.label = [mzLabel,' / ', mzUnit];
clustersOut.plot.z.value = codeSize;
        
%% NESTED FUNCTIONS
end

%% SUB FUNCTIONS
% 1. INITFUNCTION
% Function that get the input argument and check for errors
% DIL errors check should be used
function [parameters, options] = ...
    initFunction(narginIn, finneeStc, address, HL, minInt, vararginIn)

% 1.1. Defined default parameters
options.display = 1;
parameters.minReptMotif = 2;

% 1.2. Check obligatory parameters
if narginIn <  4 % Check the number of input parameters
    error('myApp:argChk', ...
        ['Wrong number of input arguments. \n', ...
        'Type help plotTrace for more information']);
elseif ~ischar(address)
    error('myApp:argChk', ...
        ['DATASET shoud be a string. \n', ...
        'Type help MSdata2struct for more information']);
elseif ~isstruct(finneeStc)
    error('myApp:argChk', ...
        ['finneeStc shoud be a structure. \n', ...
        'Type help MSdata2struct for more information'])
end
parameters.minIntensity = minInt;
parameters.HL = HL;

% 1.3. Decifer address
list = strsplit(address, '@');
tgtDataset = str2double(list{2});
tgtHStr =  str2double(list{1});
parameters.dataset = tgtDataset; 
parameters.HStruct = tgtHStr;

% 1.4. Check  optional parameters
if  narginIn > 4
    SFi = 1;
    length(vararginIn)
    while SFi <= length(vararginIn)
        switch vararginIn{SFi}
            case 'minPIPperCluster'
                parameters.minReptMotif = vararginIn{SFi+1};
                SFi = SFi + 2;
            case 'displayOff'
                options.display = 0;
            otherwise
                error('myApp:argChk', ...
                    [vararginIn{SFi} ' is not a recognized PropertyName'])
        end
    end
end
end

