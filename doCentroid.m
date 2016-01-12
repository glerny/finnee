function finneeStc = doCentroid (finneeStc, dataset, varargin)
%% DESCRIPTION
% 1. INTRODUCTION
% DOCENTROID is a function that is used to calculate the the centroid
% spectrum dataset from a profile spectrum dataset. With this function the
% centroid spectrum are calculated using local maxima.
%
% 2. PARAMETERS:
%   .required. DOCENTROID requires at least 2 parameters
%       finneeStc
%           is the finnee structure that contain information about the run
%           and link and indexation of the associated dat file. The
%           strcuture should have been create by function such as
%           MZML2STRUCT
%       dataset
%           datasetis the indice to the targeted dataset (i.e. in
%           finneeStc.dataset{m}, where m is the target dataset). The
%           dataset can be a 'centroid spectrum' a 'profile spectrum' or a
%           'ionic profile' dataset
%
%   .optionals. VARARGIN describes the optional paramters.
%       'Method' followed my the name of the centroid method
%           gle1' : default method 
%                   find limits via local minima and maxima
%                   the peak picking algorithms and output can be fund at
%                   /subfunction/peakPickingMS_gle1.m
%       'Tmin:Tmax'
%           Not implemented yet
%       'MZmin:MZmax'
%           Not implemented yet
%
% 3. EXAMPLES:
%
%
% 4. COPYRIGHT
% Copyright 2014-2015 G. Erny (guillaume@fe.up.pt), FEUP, Porto, Portugal
%

%% CORE OF THE FUNCTION
% 1. INITIALISATION
disp('hello')
info.functionName = 'doCentroid';
info.description{1} = 'do a ''centroid spectroid'' dataset from a  ''profile spectrum'' dataset';
info.matlabVersion = '8.5.0.197613 (R2015a)';
info.version = '13/07/2015_gle01';
info.ownerContact = 'guillaume@fe.up.pt';
[parameters, options] = ...
    initFunction(nargin,  finneeStc, dataset, varargin );
%INITFUNCTION used to verify the entries and load the optional and
% compulsory parameters

m = parameters.dataset;
finneeStc.dataset{end+1} = finneeStc.dataset{m};
if isfield(finneeStc.dataset{end}, 'trace')
    finneeStc.dataset{end} = rmfield(finneeStc.dataset{end}, 'trace');
end
finneeStc.dataset{end}.infoFunctionUsed.parameters.decimals = ...
    parameters.decimals;

[mzMin, intMin] = deal(inf); [mzMax, intMax] = deal(0);
axeX = []; TICP = [];  BPP = []; mzBPP = []; MSIndex = [];
datasetType = 'centroid spectrum';

fidReadDat = fopen(finneeStc.dataset{m}.description.path2DatFile, 'a+b');
% 2. CHECKING THE DATA TYPE
switch finneeStc.dataset{m}.description.dataFormat
    case 'profile spectrum'
        % 3. GETTING EACH SCAN
        
        index = finneeStc.dataset{m}.description.axe;
        fseek(fidReadDat, index(1), 'bof');
        axeX = fread(fidReadDat, [(index(2)-index(1))/(index(3)*8), index(3)], 'double');
        indTimeStt = findCloser(parameters.xMin, axeX);
        indTimeEnd = findCloser(parameters.xMax, axeX);
        
        for ii = indTimeStt:indTimeEnd
            disp(['processing scan ', num2str(ii), ' out of ',...
                num2str(length(axeX)), ' scans'])
            index = finneeStc.dataset{m}.description.index2DotDat(ii, :);
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
        warning('on')
        
        save2struc()
        fclose(fidReadDat);
    otherwise
        fclose(fidReadDat);
        error('dataset should be ''profile spectrum'' dataset')
end

save(fullfile(finneeStc.infoFunctionUsed.parameters.folderOut, ...
    [finneeStc.infoFunctionUsed.parameters.fileID '.fin']), 'finneeStc', '-mat')


%% NESTED FUNCTIONS
% 1. SAVE2STRUCT
% save results to the finnee structure
    function save2struc()
        finneeStc.dataset{end}.infoFunctionUsed.info = info;
        finneeStc.dataset{end}.infoFunctionUsed.parameters = parameters;
        
        % 1.1. Save description of dataset
        finneeStc.dataset{end}.description.mzStart = mzMin;
        finneeStc.dataset{end}.description.mzEnd = mzMax;
        finneeStc.dataset{end}.description.intMin = intMin;
        finneeStc.dataset{end}.description.intMax = intMax;
        finneeStc.dataset{end}.description.timeStart = axeX(indTimeStt);
        finneeStc.dataset{end}.description.timeEnd = axeX(indTimeEnd);
        finneeStc.dataset{end}.description.dataFormat = datasetType;
        finneeStc.dataset{end}.description.index2DotDat = MSIndex;
        
        % Record axeX in the *dat* file the rest in the *tra* file
        fseek(fidReadDat, 0, 'eof');
        finneeStc.dataset{end}.description.axe = [ftell(fidReadDat), 0, 1];
        axeX = axeX(indTimeStt:indTimeEnd);
        fwrite(fidReadDat, axeX, 'double');
        finneeStc.dataset{end}.description.axe(2) = ftell(fidReadDat);
        
        % 1.2. Record profiles
        fseek(fidReadDat, 0, 'eof');
        
        % ** TICP
        finneeStc.dataset{end}.trace{1}.infoFunctionUsed.info = info;
        finneeStc.dataset{end}.trace{1}.infoFunctionUsed.parameters = parameters;
        finneeStc.dataset{end}.trace{1}.description.name = ...
            ['Total Ion Current Profile (dataset ', ...
            num2str(length(finneeStc.dataset)), ')'];
        finneeStc.dataset{end}.trace{1}.description.dateOfCreation = clock;
        finneeStc.dataset{end}.trace{1}.description.plotType = 'profile';
        finneeStc.dataset{end}.trace{1}.description.axeX.label = ...
            finneeStc.dataset{end}.description.timeLabel;
        finneeStc.dataset{end}.trace{1}.description.axeX.unit = ...
            finneeStc.dataset{end}.description.timeUnit;
        finneeStc.dataset{end}.trace{1}.description.axeY.label = ...
            finneeStc.dataset{end}.description.intLabel;
        finneeStc.dataset{end}.trace{1}.description.axeY.unit = ...
            finneeStc.dataset{end}.description.intUnit;
        finneeStc.dataset{end}.trace{1}.index2DotDat = ...
            [ftell(fidReadDat), 0, 2];
        fwrite(fidReadDat, [axeX TICP'], 'double');
        finneeStc.dataset{end}.trace{1}.index2DotDat(2) = ftell(fidReadDat);
        
        % ** BPP
        finneeStc.dataset{end}.trace{2}.infoFunctionUsed.info = info;
        finneeStc.dataset{end}.trace{2}.infoFunctionUsed.parameters = parameters;
        finneeStc.dataset{end}.trace{2}.description.name = ...
            ['Base Peak Profile (dataset ', ...
            num2str(length(finneeStc.dataset)), ')'];
        finneeStc.dataset{end}.trace{2}.description.dateOfCreation = clock;
        finneeStc.dataset{end}.trace{2}.description.plotType = 'profile';
        finneeStc.dataset{end}.trace{2}.description.axeX.label = ...
            finneeStc.dataset{end}.description.timeLabel;
        finneeStc.dataset{end}.trace{2}.description.axeX.unit = ...
            finneeStc.dataset{end}.description.timeUnit;
        finneeStc.dataset{end}.trace{2}.description.axeY.label = ...
            finneeStc.dataset{end}.description.intLabel;
        finneeStc.dataset{end}.trace{2}.description.axeY.unit = ...
            finneeStc.dataset{end}.description.intUnit;
        finneeStc.dataset{end}.trace{2}.index2DotDat = ...
            [ftell(fidReadDat), 0, 2];
        fwrite(fidReadDat, [axeX BPP'], 'double');
        finneeStc.dataset{end}.trace{2}.index2DotDat(2) = ftell(fidReadDat);
        
        % ** mzBPP
        finneeStc.dataset{end}.trace{3}.infoFunctionUsed.info = info;
        finneeStc.dataset{end}.trace{3}.infoFunctionUsed.parameters = parameters;
        finneeStc.dataset{end}.trace{3}.description.name = ...
            ['m/z @ Base Peak (dataset ', ...
            num2str(length(finneeStc.dataset)), ')'];
        finneeStc.dataset{end}.trace{3}.description.dateOfCreation = clock;
        finneeStc.dataset{end}.trace{3}.description.plotType = 'profile';
        finneeStc.dataset{end}.trace{3}.description.axeX.label = ...
            finneeStc.dataset{end}.description.timeLabel;
        finneeStc.dataset{end}.trace{3}.description.axeX.unit = ...
            finneeStc.dataset{end}.description.timeUnit;
        finneeStc.dataset{end}.trace{3}.description.axeY.label = ...
            finneeStc.dataset{end}.description.intLabel;
        finneeStc.dataset{end}.trace{3}.description.axeY.unit = ...
            finneeStc.dataset{end}.description.intUnit;
        finneeStc.dataset{end}.trace{3}.index2DotDat = ...
            [ftell(fidReadDat), 0, 2];
        fwrite(fidReadDat, [axeX mzBPP'], 'double');
        finneeStc.dataset{end}.trace{3}.index2DotDat(2) = ftell(fidReadDat);
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
parameters.decimals = 5;
parameters.method = 'gle1';

% 2. Check for option
if  narginIn > 2
    SFi = 1;
    length(vararginIn)
    while SFi <= length(vararginIn)
        switch vararginIn{SFi}
            case 'Method'
                parameters.method = vararginIn{SFi+1};
                SFi = SFi +2;
            otherwise
                error('myApp:argChk', ...
                    [vararginIn{SFi} ' is not a recognized PropertyName'])
        end
    end
end
end
