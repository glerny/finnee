function dataOut = doClustersPlot(finneeStc, address, HL, minInt, varargin)
%% DESCRIPTION
% 1. INTRODUCTION
% DOCLUSTERPLOT allows to select clusters at given hierarchical level
% defined as the minimal correlation coefficient between pure ion profiles
% (PIP) within each cluster.
%
% 2. PARAMETERS:
%   .required. DOCLUSTERPLOT requires at least 4 parameters
%       finneeStc
%           is the finnee structure that contain information about the run
%           and link and indexation of the associated dat file. The
%           strcuture should have been create by function such as 
%           MZML2STRUCT
%       addresss
%           defined the target cluster analysis (HACA) within the target 
%           dataset. The format should be 'HACA@dataset'.
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

%   .optionals. VARARGIN describes the optional paramters.  
%       'displayOff' 
%           Allow to avoid any display in the matlab command window
%       'minPIPperCluster' followed by an integer (default 2)
%           Only record clusters with 'minPIPperCluster' PIP
%
% 3. EXAMPLES:
%	dataOut = doClustersPlot(finneeStc, '1@2', 0.95, 500)
%
% 4. COPYRIGHT
%   Copyright 2014-2015 G. Erny (guillaume@fe.up.pt), FEUP, Porto, Portugal

%% CORE OF THE FUNCTION
% 1. INITIALISATION
info.functionName = 'doHACA';
info.description{1} = 'Calculate the herarchical structure of '' dataset';
info.matlabVersion = '8.5.0.197613 (R2015a)';
info.version = '10/07/2015_gle01';
info.ownerContact = 'guillaume@fe.up.pt';

[parameters, options] = ...
    initFunction(nargin, finneeStc, address, HL, minInt, varargin);
%INITFUNCTION - sub function used to verify the entries and load the optional and
% compulsory and optional parameters

m = parameters.dataset;
n = parameters.HStruct;

fidReadDat = fopen(finneeStc.dataset{m}.description.path2DatFile, 'rb');
fom.label = finneeStc.dataset{m}.description.fom.label;
index = finneeStc.dataset{m}.description.fom.data;
fseek(fidReadDat, index(1), 'bof');
fom.data = fread(fidReadDat, ...
    [(index(2)-index(1))/(index(3)*8), index(3)], 'double');
index = finneeStc.dataset{m}.description.axe;
fseek(fidReadDat, index(1), 'bof');
axeX = fread(fidReadDat, [(index(2)-index(1))/(index(3)*8), ...
    index(3)], 'double');
        
formatSpec = ['%.',...
    num2str(finneeStc.dataset{m}.infoFunctionUsed.parameters.decimals), 'f'];
% TOBECHANGED

clusters.step = ...
    finneeStc.dataset{m}.clusterAnalysis{n}.HStructure.step;

% 2. FIND THE HL AND LOAD TARGETED CLUSTERS
targetStep = [];
for ii = 1:length(clusters.step)
    if clusters.step{ii}.HL < parameters.HL && isempty(targetStep)
        targetStep = ii-1;
    end
end
if isempty(targetStep), targetStep = ii ; end


index = finneeStc.dataset{m}.clusterAnalysis{n}...
    .HStructure.step{targetStep}.index2DotDat;
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

figure
PIPinClust = 0;
for ii = 1:length(selClusters)
    % Calculate SumSIPinCl
    sumPIPinCl = zeros(length(axeX), 1);
    for jj = 1:length(selClusters{ii})
        index = finneeStc.dataset{m}.description...
            .index2DotDat(selClusters{ii}(jj), :)
        fseek(fidReadDat, index(1), 'bof');
        PIP = fread(fidReadDat, [(index(2)-index(1))/(index(3)*8), ...
            index(3)], 'double');
        sumPIPinCl(PIP(:,1)) = sumPIPinCl(PIP(:,1)) + PIP(:,2);
    end
    if max(fom.data(selClusters{ii}, 5)) >= parameters.minIntensity
        PIPinClust = PIPinClust +  length(selClusters{ii});
        sumCluSIP = sumCluSIP + sumPIPinCl;
        sumPIP = sumPIP + sumPIPinCl;
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
        subplot(2,2,2)
        hold on
        plot(axeX, sumPIPinCl, 'k')
        subplot(2,2,4)
        hold on
        plot(axeX,sumPIPinCl/max(sumPIPinCl), 'k')
    else
        sumPIP = sumPIP + sumPIPinCl;
    end
end
fclose(fidReadDat);

subplot(2,2,1)
plot(axeX, sumCluSIP, 'k')
title('Sum of single ion profiles');
xlabel([finneeStc.dataset{m}.description.timeLabel,' / ', ...
    finneeStc.dataset{m}.description.timeUnit]);
ylabel([finneeStc.dataset{m}.description.intLabel,' / ', ...
    finneeStc.dataset{m}.description.intUnit]);

subplot(2,2,2)
title('Most intensed PIP in each cluster');
xlabel([finneeStc.dataset{m}.description.timeLabel,' / ', ...
    finneeStc.dataset{m}.description.timeUnit]);
ylabel([finneeStc.dataset{m}.description.intLabel,' / ', ...
    finneeStc.dataset{m}.description.intUnit]);
hold off
subplot(2,2,4)
title('Most intensed PIP in each cluster (normalised intensities)');
xlabel([finneeStc.dataset{m}.description.timeLabel,' / ', ...
    finneeStc.dataset{m}.description.timeUnit]);
ylabel([finneeStc.dataset{m}.description.intLabel,' / ', ...
    finneeStc.dataset{m}.description.intUnit]);
hold off
subplot(2,2,3)
title('Residual all PIP - PIP in clusters');
plot(axeX, sumPIP-sumCluSIP, 'r')
xlabel([finneeStc.dataset{m}.description.timeLabel,' / ', ...
    finneeStc.dataset{m}.description.timeUnit]);
ylabel([finneeStc.dataset{m}.description.intLabel,' / ', ...
    finneeStc.dataset{m}.description.intUnit]);hold off

% 4. CLUSTERS PLOT
listMZ2plot3D = fomClusters.data(:,6);
listI2plot3D =  fomClusters.data(:,9);
listT2plot3D = fomClusters.data(:,2);
figure('units','normalized','outerposition',[0 0 1 1], 'Name', ... % this does make it full size on my computer don't knwo why
    ['3D display of HACBotUp of ', finneeStc.infoFunctionUsed.parameters.fileID, ' @ ',...
    address, ' with deltaM1 <= ', num2str(parameters.HL), ...
    ' repetition of motifs >= ' , num2str(parameters.minReptMotif), ...
    ' and intenisty of clusters >= ', num2str(parameters.minIntensity)])

subplot(5, 1, 1) % TICP using only pure ions profiles in clusters
plot(axeX, sumCluSIP, 'k')
title('Sum of single ion profiles');
xlabel([finneeStc.dataset{m}.description.timeLabel,' / ', ...
    finneeStc.dataset{m}.description.timeUnit]);
ylabel([finneeStc.dataset{m}.description.intLabel,' / ', ...
    finneeStc.dataset{m}.description.intUnit]);
v1 = axis;

subplot(5,1, 2:5) % CLUSTERS PLOT
codeSize = int32(listI2plot3D/max(listI2plot3D)*(500) + 5);
hold on
scatter(listT2plot3D, listMZ2plot3D, ...
    codeSize, 'k');
assignin('base', 'listT2plot3D', listT2plot3D)
assignin('base', 'listMZ2plot3D', listMZ2plot3D)
assignin('base', 'codeSize', codeSize)
title('Clusters Plot');
xlabel([finneeStc.dataset{m}.description.timeLabel,' / ', ...
    finneeStc.dataset{m}.description.timeUnit]);
ylabel([finneeStc.dataset{m}.description.mzLabel,' / ', ...
    finneeStc.dataset{m}.description.mzUnit]);
v2 = axis;
axis([v1(1) v1(2) v2(3) v2(4)])

% DataOut is used by PLOTCLUSTER to analyse each cluster individually
dataOut.results.path2DatFile = ...
    finneeStc.dataset{m}.description.path2DatFile;
dataOut.results.index2DotDat = ...
    finneeStc.dataset{m }.description.index2DotDat;
dataOut.results.axeX = axeX;
dataOut.results.fom = fom;
dataOut.results.fomClusters = fomClusters;
dataOut.results.selClusters = selClusters;
dataOut.results.formatSpec = formatSpec;
dataOut.info = info;
dataOut.parameters = parameters;
        
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

