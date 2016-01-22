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
%       'noFig' 
%           No figures displayed
%       'indice'  
%           Work only with 'profile spectrum' dataset. If indice is used,
%           massInt is the indice in the mass axe.
%
% 3. EXAMPLES:
%       spectraOut = getProfile(finneeStc, 1, [500.23 500.25])
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
fmt = finneeStc.info.parameters.prec4mz;
        
% 1. Checking the data type
switch finneeStc.dataset{m}.description.dataFormat
    case 'profile spectrum'
        
        %2.1. FInd limits
        index = finneeStc.dataset{m}.indexInDat(1, :);
        fseek(fidReadDat, index(1), 'bof');
        MS = ...
            fread(fidReadDat, [(index(2)-index(1))/(index(3)*8), index(3)], 'double');
        if parameters.mzMin ~=  parameters.mzMax
            if options.indice
                indMzStt = parameters.mzMin;
                indMzEnd = parameters.mzMax;
            else
                indMzStt = findCloser(parameters.mzMin, MS(:,1));
                indMzEnd = findCloser(parameters.mzMax,  MS(:,1));
            end
            profileOut.title = ['Extracted ion profile ( m/z = ', ...
                num2str(MS(indMzStt, 1), fmt),':', ...
                num2str(MS(indMzEnd, 1), fmt), '); (',...
                'dataset ',  num2str(m), ')'];
        else
            if options.indice
                ind = parameters.mzMin;
            else
                ind = findCloser( parameters.mzMin, MS(:,1));
            end
            profileOut.title = ['Extracted ion profile ( m/z = ', ...
                num2str(MS(ind, 1),fmt),'); (',...
                'dataset ',  num2str(m), ')'];
        end
        profileOut.plotType = 'profile';
        
        profileOut.data = axeX;
        profileOut.data(:,2) = 0;
        %2.2. getting each MS spectra and caulcaulating the profile
        for ii = 1:length(axeX)
            index = finneeStc.dataset{m}.indexInDat(ii, :);
            fseek(fidReadDat, index(1), 'bof');
            MS = ...
                fread(fidReadDat, [(index(2)-index(1))/(index(3)*8), index(3)], 'double');
            if parameters.mzMin ~=  parameters.mzMax
                MS = MS(indMzStt:indMzEnd, :);
            else
                MS = MS(ind, :);
            end
            if isempty(MS)
                profileOut.data(ii, 2) = 0;
            else
                profileOut.data(ii, 2) = sum(MS(:,2));
            end
        end
        
    case 'centroid spectrum'
        if parameters.mzMin ~=  parameters.mzMax
            profileOut.title = ['Extracted ion profile ( m/z = ', ...
                num2str( parameters.mzMin, fmt),':', ...
                num2str(parameters.mzMax, fmt), '); (',...
                'dataset ',  num2str(m), ')'];
        else
            profileOut.title = ['Extracted ion profile ( m/z = ', ...
                num2str(parameters.mzMin),'); (',...
                'dataset ',  num2str(m), ')'];
        end
        profileOut.plotType = 'profile';
        profileOut.data = axeX;
        profileOut.data(:,2) = 0;
        
        %2.2. getting each MS spectra and caulcaulating the profile
        for ii = 1:length(axeX)
            index = finneeStc.dataset{m}.indexInDat(ii, :);
            fseek(fidReadDat, index(1), 'bof');
            MS = ...
                fread(fidReadDat, [(index(2)-index(1))/(index(3)*8), index(3)], 'double');
            ind2rem = MS(:,1) < parameters.mzMin |...
                MS(:,1) > parameters.mzMax;
            MS(ind2rem, :) = [];
            if isempty(MS)
                profileOut.data(ii, 2) = 0;
            else
                profileOut.data(ii, 2) = sum(MS(:,2));
            end
        end
    case 'ionic profile'
         profileOut.plotType = 'profile';
         if parameters.mzMin ~=  parameters.mzMax
            profileOut.title = ['Extracted ion profile ( m/z = ', ...
                num2str( parameters.mzMin, fmt),':', ...
                num2str(parameters.mzMax, fmt), '); (',...
                'dataset ',  num2str(m), ')'];
        else
            profileOut.title = ['Extracted ion profile ( m/z = ', ...
                num2str(parameters.mzMin),'); (',...
                'dataset ',  num2str(m), ')'];
        end
        
        profileOut.data = axeX;
        profileOut.data(:,2) = 0;
        %2.3 getting each PIP
        for ii = 1:length(finneeStc.dataset{m}.indexInDat(:,1))
            index = finneeStc.dataset{m}.indexInDat(ii, :);
            fseek(fidReadDat, index(1), 'bof');
            PIP = fread(fidReadDat, [(index(2)-index(1))/(index(3)*8), ...
                index(3)], 'double');
            ind2rem = PIP(:,3) <  parameters.mzMin |...
                PIP(:,3) > parameters.mzMax;
            PIP(ind2rem, :) = [];
            if ~isempty(PIP)
                profileOut.data(PIP(:,1),2) = profileOut.data(PIP(:,1),2) + PIP(:,2);
            end
        end
    otherwise
        DIL
end
fclose(fidReadDat);

profileOut.axes.axeX.label = finneeStc.dataset{m}.axes.time.label;
profileOut.axes.axeX.unit = finneeStc.dataset{m}.axes.time.unit;
profileOut.axes.axeY.label = finneeStc.dataset{m}.axes.intensity.label;
profileOut.axes.axeY.unit = finneeStc.dataset{m}.axes.intensity.unit;

% 3. Plot if requested
if options.display
    switch profileOut.plotType
        case 'profile'
            plot(profileOut.data(:,1), profileOut.data(:,2));
            title(profileOut.title);
            xlabel([profileOut.axes.axeX.label, ' / ', profileOut.axes.axeX.unit]);
            ylabel([profileOut.axes.axeY.label, ' / ', profileOut.axes.axeY.unit]);
        case 'stem'
            stem(profileOut.data(:,1), profileOut.data(:,2), 'Marker', 'none');
            title(profileOut.title);
            xlabel([profileOut.axes.axeX.label, ' / ', profileOut.axes.axeX.unit]);
            ylabel([profileOut.axes.axeY.label, ' / ', profileOut.axes.axeY.unit]);
        case 'bar'
            bar(profileOut.data(:,1), profileOut.data(:,2));
            title(profileOut.title);
            xlabel([profileOut.axes.axeX.label, ' / ', profileOut.axes.axeX.unit]);
            ylabel([profileOut.axes.axeY.label, ' / ', profileOut.axes.axeY.unit]);
    end
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
            case 'noFig'
                options.display = 0;
                SFi = SFi + 1;
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
