function finneeStc = doCentroid (finneeStc, dataset, varargin)
%% DESCRIPTION
% 1. INTRODUCTION
% DOCENTROID is a function that is used to calculate the the centroid
% spectrum dataset from a profile spectrum dataset. With this function the
% centroid spectrum are calculated using local maxima.
%
% 2. INPUT PARAMETERS
%   .required. DOCENTROID requires at least 2 parameters
%       finneeStc
%           is the finnee structure that contain information about the run
%           and link and indexation of the associated dat file. 
%       dataset
%           dataset is the indice to the targeted dataset (i.e. in
%           finneeStc.dataset{m}, where m is the target dataset). The
%           dataset sould be in a 'profile spectrum' format
%
%   .optionals. VARARGIN describes the optional paramters.
%       'method' followed my the name of the centroid method
%           gle1' : default method 
%                   find  maxima. The peak picking algorithms and output
%                   can be fund at peakPickingMS_gle1.m
%
% 3. OUTPUT PARAMETERS
%   .finneeStc 
%       Target Finnee structure 
%
% 4. EXAMPLES
%    finneeStc = doCentroid(finneeStc);
%
% 5. COPYRIGHT
% Copyright 2015-2016 G. Erny (guillaume@fe.up.pt), FEUP, Porto, Portugal
%

%% CORE OF THE FUNCTION
% 1. INITIALISATION
info.function.functionName = 'doCentroid';
info.function.description{1} = 'do a ''centroid spectrum'' dataset from a  ''profile spectrum'' dataset';
info.function.matlabVersion = '8.5.0.197613 (R2015a)';
info.function.version = '14/01/2016';
info.function.ownerContact = 'guillaume@fe.up.pt';
[parameters, options] = ...
    initFunction(nargin,  finneeStc, dataset, varargin );
%INITFUNCTION used to verify the entries and load the optional and
% compulsory parameters

m = parameters.dataset;

[mzMin, intMin] = deal(inf); [mzMax, intMax] = deal(0);
TICP = [];  BPP = []; mzBPP = []; MSIndex = [];
datasetType = 'centroid spectrum';

fidReadDat = fopen(finneeStc.path2dat, 'a+b');
% 2. CHECKING THE DATA TYPE
switch finneeStc.dataset{m}.description.dataFormat
    case 'profile spectrum'
        % 3. GETTING EACH SCAN
        
        axeX = finneeStc.dataset{m}.axes.time.values;
        indTimeStt = findCloser(parameters.xMin, axeX);
        indTimeEnd = findCloser(parameters.xMax, axeX);
        
        for ii = indTimeStt:indTimeEnd
            disp(['processing scan ', num2str(ii), ' out of ',...
                num2str(length(axeX)), ' scans'])
            index = finneeStc.dataset{m}.indexInDat(ii, :);
            fseek(fidReadDat, index(1), 'bof');
            MS = ...
                fread(fidReadDat, [(index(2)-index(1))/(index(3)*8), index(3)], 'double');
            ind2rem = MS(:,1) < parameters.mzMin |...
                MS(:,1) > parameters.mzMax;
            MS(ind2rem, :) = [];
            if ~isempty(MS)
                
                % 4. DOING THE PEAKPICKING
                switch parameters.method
                    case 'gle1'
                        ppMS = peakPickingMS_gle1(MS);
                        
                        % 5. CALCULATE PROFILE AND SAVE DATA
                        % !NOTE! the code is set for MZ at peak apex 
                        
                        if mzMin > min(ppMS(:,1)), mzMin = min(ppMS(:,1)); end
                        if mzMax < max(ppMS(:,1)), mzMax = max(ppMS(:,1)); end
                        if intMin > min(ppMS(:,2)), intMin = min(ppMS(:,2)); end
                        if intMax < max(ppMS(:,2)), intMax = max(ppMS(:,2)); end
                        
                        if isempty(ppMS)
                            TICP(ii) = 0;
                            BPP(ii) = 0;
                            mzBPP(ii) = 0;
                            fseek(fidReadDat, 0, 'eof')
                            MSIndex(ii, :) = [ftell(fidReadDat), ftell(fidReadDat), 2];
                        else
                            TICP(ii) = sum(ppMS(:,2));
                            [BPP(ii), indMax] = max(ppMS(:,2));
                            mzBPP(ii) = ppMS(indMax,1);
                            fseek(fidReadDat, 0, 'eof');
                            MSIndex(ii, :) = [ftell(fidReadDat), 0, 2];
                            fwrite(fidReadDat, ppMS, 'double');
                            MSIndex(ii, 2) = ftell(fidReadDat);
                        end
                    otherwise
                        error('This method has not been implemented in this version')
                end
                
               
            end
        end
        save2struc()
        fclose(fidReadDat);
    otherwise
        fclose(fidReadDat);
        error('dataset should be ''profile spectrum'' dataset')
end

save(fullfile(finneeStc.info.parameters.folderOut, ...
    [finneeStc.info.parameters.fileID '.mat']), 'finneeStc')


%% NESTED FUNCTIONS
% 1. SAVE2STRUCT
% save results to the finnee structure
    function save2struc()
        finneeStc.dataset{end+1}.name = ...
            ['''centroid spectrum'' dataset of the ''profile spectrum'' dataset ', ...
            num2str(m)];
        finneeStc.dataset{end}.dateOfCreation = datetime;
        finneeStc.dataset{end}.info = info;
        finneeStc.dataset{end}.info.parameters = parameters;
        finneeStc.dataset{end}.info.errors = {};

        timeLabel = finneeStc.dataset{m}.axes.time.label;
        timeUnit = finneeStc.dataset{m}.axes.time.unit;
        mzLabel = finneeStc.dataset{m}.axes.mz.label;
        mzUnit = finneeStc.dataset{m}.axes.mz.unit;
        intLabel = finneeStc.dataset{m}.axes.intensity.label;
        intUnit = finneeStc.dataset{m}.axes.intensity.unit;
        
        % 1.1. Save description of dataset
        finneeStc.dataset{end}.description.mzStart = mzMin;
        finneeStc.dataset{end}.description.mzEnd = mzMax;
        finneeStc.dataset{end}.description.intMin = intMin;
        finneeStc.dataset{end}.description.intMax = intMax;
        finneeStc.dataset{end}.description.timeStart = axeX(indTimeStt);
        finneeStc.dataset{end}.description.timeEnd = axeX(indTimeEnd);
        finneeStc.dataset{end}.description.dataFormat = datasetType;
        finneeStc.dataset{end}.indexInDat = MSIndex;
        
        finneeStc.dataset{end}.axes.time.values = axeX;
        finneeStc.dataset{end}.axes.time.label = timeLabel;
        finneeStc.dataset{end}.axes.time.unit = timeUnit;
        finneeStc.dataset{end}.axes.mz.values = [];
        finneeStc.dataset{end}.axes.mz.label = mzLabel;
        finneeStc.dataset{end}.axes.mz.unit = mzUnit;
        finneeStc.dataset{end}.axes.intensity.values = [];
        finneeStc.dataset{end}.axes.intensity.label = intLabel;
        finneeStc.dataset{end}.axes.intensity.unit = intUnit;
        
        % 1.2. Record profiles
        % ** TICP
        finneeStc.dataset{end}.trace{1}.name = ...
            ['Total Ion Current Profile (dataset ', ...
            num2str(length(finneeStc.dataset)), ')'];
        finneeStc.dataset{end}.trace{1}.dateOfCreation = datetime;
        finneeStc.dataset{end}.trace{1}.code = 'TIP';
        finneeStc.dataset{end}.trace{1}.plotType = 'profile';
        finneeStc.dataset{end}.trace{1}.axeX.label = timeLabel;
        finneeStc.dataset{end}.trace{1}.axeX.unit = timeUnit;
        finneeStc.dataset{end}.trace{1}.axeY.label = intLabel;
        finneeStc.dataset{end}.trace{1}.axeY.unit = intUnit;
        finneeStc.dataset{end}.trace{1}.indexInDat  = [ftell(fidReadDat), 0, 2];
        fwrite(fidReadDat, [axeX TICP'], 'double');
        finneeStc.dataset{end}.trace{1}.indexInDat(2) = ftell(fidReadDat);
        
        % ** BPP
        finneeStc.dataset{end}.trace{2}.name = ...
            ['Base Peak Profile (dataset ', ...
            num2str(length(finneeStc.dataset)), ')'];
        finneeStc.dataset{end}.trace{2}.dateOfCreation = datetime;
        finneeStc.dataset{end}.trace{1}.code = 'BPP';
        finneeStc.dataset{end}.trace{2}.plotType = 'profile';
        finneeStc.dataset{end}.trace{2}.axeX.label = timeLabel;
        finneeStc.dataset{end}.trace{2}.axeX.unit = timeUnit;
        finneeStc.dataset{end}.trace{2}.axeY.label = intLabel;
        finneeStc.dataset{end}.trace{2}.axeY.unit = intUnit;
        finneeStc.dataset{end}.trace{2}.indexInDat  = [ftell(fidReadDat), 0, 2];
        fwrite(fidReadDat, [axeX BPP'], 'double');
        finneeStc.dataset{end}.trace{2}.indexInDat(2) = ftell(fidReadDat);
        
        % ** mzBPP
        finneeStc.dataset{end}.trace{3}.name = ...
            ['m/z @ Base Peak (dataset ', ...
            num2str(length(finneeStc.dataset)), ')'];
        finneeStc.dataset{end}.trace{3}.dateOfCreation = datetime;
        finneeStc.dataset{end}.trace{1}.code = 'mzBPP';
        finneeStc.dataset{end}.trace{3}.plotType = 'profile';
        finneeStc.dataset{end}.trace{3}.axeX.label = timeLabel;
        finneeStc.dataset{end}.trace{3}.axeX.unit = timeUnit;
        finneeStc.dataset{end}.trace{3}.axeY.label = mzLabel;
        finneeStc.dataset{end}.trace{3}.axeY.unit = mzUnit;
        finneeStc.dataset{end}.trace{3}.indexInDat  = [ftell(fidReadDat), 0, 2];
        fwrite(fidReadDat, [axeX mzBPP'], 'double');
        finneeStc.dataset{end}.trace{3}.indexInDat(2) = ftell(fidReadDat);
    end



end


%% SUB FUNCTIONS
% 1. INITFUNCTION
% Function that get the input argument and check for errors
function [parameters, options] = ...
    initFunction(narginIn, finneeStc, dataset, vararginIn )

% 1. Check for obligatory parameters
if narginIn < 2 % check the number of input parameters
    error('myApp:argChk', ...
        ['Wrong number of input arguments. \n', ...
        'Type help plotTrace for more information']);
elseif ~isnumeric(dataset)
    error('myApp:argChk', ...
        ['DATASET shoud be a number. \n', ...
        'Type help MSdata2struct for more information']);
elseif ~isstruct(finneeStc)
    error('myApp:argChk', ...
        ['finneeStc shoud be a structure. \n', ...
        'Type help MSdata2struct for more information']);
end

parameters.mzMin = 0;
parameters.mzMax = inf;
parameters.xMin = 0;
parameters.xMax = inf;
options.text = 1;
options.save = 1;
parameters.dataset = dataset;
parameters.method = 'gle1';

% 2. Check for option
if  narginIn > 2
    SFi = 1;
    length(vararginIn)
    while SFi <= length(vararginIn)
        switch vararginIn{SFi}
            case 'method'
                parameters.method = vararginIn{SFi+1};
                SFi = SFi +2;
            otherwise
                error('myApp:argChk', ...
                    [vararginIn{SFi} ' is not a recognized PropertyName'])
        end
    end
end
end
