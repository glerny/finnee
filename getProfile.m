function [profileOut, finneeStc] = getProfile(finneeStc, dataset, massInt, varargin)
%% DESCRIPTION
% 1. INTRODUCTION
% GETPROFILE is used to generate a profile using a given m/z interval. 
% This function will update the finneeStc in case the
% profile is to be save (optional. profileOut is a 2xm array that contain
% the profile. This function work with profile or centroid spectrum, and 
% extracted ions dataset. 
%
% 2. PARAMETERS:
%   .required. GETPROFILE requires at least 3 parameters
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
%       massInt
%           Define the m/z interval to be used. It can be a single value
%           (closest m/z value) or a 2x1 array [mzMin mzMax]
%
%   .optionals. VARARGIN describes the optional paramters.  
%       'timeInt' follow by a 2x1 array
%           Allow to only take the data points with m/z are within the
%           defined interval
%       'noFig' 
%           No figures displayed
%       'indice'  
%           Work only with 'profile spectrum' dataset. If indice is used,
%           massInt is the indice in the mass axe.
%
% 3. EXAMPLES:
%       spectraOut = getProfile(finneeStc, 1, [500.23 500.25], 'timeInt', {5 10])
%
% 4. COPYRIGHT
% Copyright 2015-2016 G. Erny (guillaume@fe.up,pt), FEUP, Porto, Portugal

%% CORE OF THE FUNCTION
% 1. INITIALISATION
info.function.functionName =  'getProfile';
info.function.description{1} =  'get the profile between a m/z interval';
info.function.matlabVersion = '8.5.0.197613 (R2015a)';
info.function.version = '14/01/2016';
info.function.ownerContact = 'guillaume@fe.up.pt';
[parameters, options] = ...
    initFunction(nargin, finneeStc, dataset, massInt, varargin );
%INITFUNCTION used to verify the entries and load the optional and
% complusory parameters

m = parameters.dataset;
fidReadDat = fopen(finneeStc.path2dat, 'rb');

%% FUNCTION CORE

axeX = finneeStc.dataset{m}.axes.time.values;
profileOut.data = axeX;
profileOut.data(:,2) = 0;
        
% 1. Checking the data type
switch finneeStc.dataset{m}.description.dataFormat
    case 'profile spectrum'
        plotType = 'profile';
        
        %2.1. getting each MS spectra and caulcaulating the profile
        for ii = 1:length(axeX)
            index = finneeStc.dataset{m}.description.index2DotDat(ii, :);
            fseek(fidReadDat, index(1), 'bof');
            MS = ...
                fread(fidReadDat, [(index(2)-index(1))/(index(3)*8), index(3)], 'double');
            if parameters.mzMin ~=  parameters.mzMax
                ind2rem = MS(:,1) < parameters.mzMin |...
                    MS(:,1) > parameters.mzMax;
                MS(ind2rem, :) = [];
            else
                if options.indice
                    ind = parameters.mzMin;
                else
                    ind = findCloser( parameters.mzMin, MS(:,1));
                end
                MS = MS(ind, :);
            end
            if isempty(MS)
                profileOut(ii, 2) = 0;
            else
                switch parameters.calc
                    case 'sum'
                        profileOut(ii, 2) = sum(MS(:,2));
                    case 'max'
                        profileOut(ii, 2) = max(MS(:,2));
                    case 'freq'
                        profileOut(ii, 2) = length(MS(:,2));
                end
            end
        end
        
    case 'centroid spectrum'
        plotType = 'profile';
        
        %2.2. getting each MS spectra and caulcaulating the profile
        for ii = 1:length(axeX)
            index = finneeStc.dataset{m}.description.index2DotDat(ii, :);
            fseek(fidReadDat, index(1), 'bof');
            MS = ...
                fread(fidReadDat, [(index(2)-index(1))/(index(3)*8), index(3)], 'double');
            ind2rem = MS(:,1) < parameters.mzMin |...
                MS(:,1) > parameters.mzMax;
            MS(ind2rem, :) = [];
            if isempty(MS)
                profileOut(ii, 2) = 0;
            else
                switch parameters.calc
                    case 'sum'
                        profileOut(ii, 2) = sum(MS(:,2));
                    case 'max'
                        profileOut(ii, 2) = max(MS(:,2));
                    case 'freq'
                        profileOut(ii, 2) = length(MS(:,2));
                end
            end
        end
    case 'ionic profile'
         plotType = 'profile';
         
        %2.3 getting each PIP
        for ii = 1:length(finneeStc.dataset{m}.description.index2DotDat(:,1))
            index = finneeStc.dataset{m}.description.index2DotDat(ii, :);
            fseek(fidReadDat, index(1), 'bof');
            PIP = fread(fidReadDat, [(index(2)-index(1))/(index(3)*8), ...
                index(3)], 'double');
            ind2rem = PIP(:,3) <  parameters.mzMin |...
                PIP(:,3) > parameters.mzMax;
            PIP(ind2rem, :) = [];
            if ~isempty(PIP)
                switch parameters.calc
                    case 'sum'
                        profileOut(PIP(:,1),2) = ...
                            profileOut(PIP(:,1),2) + PIP(:,2);
                    case 'max'
                        ind2keep = PIP(:,2) > profileOut(PIP(:,1),2);
                        profileOut(PIP(ind2keep,1),2) = PIP(ind2keep,2);
                    case 'freq'
                        profileOut(PIP(:,1),2) = ...
                            profileOut(PIP(:,1),2) + 1;
                end
                profileOut(PIP(:,1),2) = profileOut(PIP(:,1),2) + PIP(:,2);
            end
        end
    otherwise
        DIL
end

if options.fid.in
else
    fclose(fidReadDat);
end


% 3. Plot and save if requested
if options.display
    strName = ['Profile with m/z between ', num2str(parameters.mzMin), ' to ', ...
        num2str(parameters.mzMax), ' (', parameters.calc, ')'];
    switch plotType
        case 'profile'
            plot(profileOut(:,1), profileOut(:,2));
        case 'stem'
            stem(profileOut(:,1), profileOut(:,2), 'Marker', 'none');
        case 'barPlot'
            bar(profileOut(:,1), profileOut(:,2), 1);
    end
    title(strName);
    xlabel([finneeStc.dataset{m}.description.timeLabel,' / ',...
        finneeStc.dataset{m}.description.timeUnit]);
    ylabel([finneeStc.dataset{m}.description.intLabel, ' / ',...
        finneeStc.dataset{m}.description.intUnit]);
end

if options.save2str
    finneeStc.dataset{m}.trace{end+1}.infoFunctionUsed.info = info;
    finneeStc.dataset{m}.trace{end}.infoFunctionUsed.parameters = parameters;
    finneeStc.dataset{m}.trace{end}.description.name = strName;
    finneeStc.dataset{m}.trace{end}.description.dateOfCreation = clock;
    finneeStc.dataset{m}.trace{end}.description.plotType = plotType;
    finneeStc.dataset{m}.trace{end}.description.axeX.label = ...
        finneeStc.dataset{m}.description.timeLabel;
    finneeStc.dataset{m}.trace{end}.description.axeX.unit = ...
        finneeStc.dataset{m}.description.timeUnit;
    finneeStc.dataset{m}.trace{end}.description.axeY.label = ...
        finneeStc.dataset{m}.description.intLabel;
    finneeStc.dataset{m}.trace{end}.description.axeY.unit = ...
        finneeStc.dataset{m}.description.intUnit;
    fidWriteTra = fopen( finneeStc.dataset{m}.description.path2DatFile, 'ab');
    fseek(fidWriteTra, 0,'eof');
    finneeStc.dataset{m}.trace{end}.index2DotDat  = ...
        [ftell(fidWriteTra), 0, 2];
    fwrite(fidWriteTra, [profileOut(:,1) profileOut(:,2)], 'double');
    finneeStc.dataset{m}.trace{end}.index2DotDat(2) = ftell(fidWriteTra);
    fclose(fidWriteTra);
    save(fullfile(finneeStc.infoFunctionUsed.parameters.folderOut, ...
        [finneeStc.infoFunctionUsed.parameters.fileID '.fin']), 'finneeStc', '-mat')
end

%% NESTED FUNCTIONS
end

%% SUB FUNCTIONS
% 1. INITFUNCTION
% Function that get the input argument and check for errors
function [parameters, options] = ...
    initFunction(narginIn, finneeStc, dataset, massInt, vararginIn )


options.display = 1;
options.save2str = 0;
options.fid.in = 0;
parameters.calc = 'sum';
options.indice = 0;

% 1.1. Check for obligatory parameters
if narginIn < 3 % check the number of input parameters
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
parameters.dataset = dataset;

if length(massInt) == 1
    [parameters.mzMin, parameters.mzMax] = deal(massInt);
elseif length(massInt) == 2
    massInt = sort(massInt);
    parameters.mzMin = massInt(1);
    parameters.mzMax = massInt(2);
else
    error('line 253')
end
parameters.xMin = 0;
parameters.xMax = inf;

% 1.2. Check for option
if  narginIn > 3
    SFi = 1;
    while SFi <= length(vararginIn)
        switch vararginIn{SFi}
            case 'timeInt'
                if length(vararginIn{SFi+1}) == 2
                    timeInt = sort(vararginIn{SFi+1});
                    parameters.xMin = timeInt(1);
                    parameters.xMax = timeInt(2);
                end
                SFi = SFi + 2;
            case 'noFig'
                options.display = 0;
                SFi = SFi + 1;
            case 'save2str'
                options.save2str = 1;
                SFi = SFi + 1;
            case 'sum'
                parameters.calc = 'sum';
                SFi = SFi + 1;
            case 'max'
                parameters.calc = 'max';
                SFi = SFi + 1;
            case 'fid'
                options.fid.in = 1 ;
                options.fid.fid = vararginIn{SFi+1};
                SFi = SFi + 2;
            case 'indice'
                options.indice = 1;
                SFi = SFi + 1;
            otherwise
                error('myApp:argChk', ...
                    [vararginIn{SFi} ' is not a recognized PropertyName'])
        end
    end
end

end