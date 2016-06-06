function spectraOut = getSpectra(finneeStc, dataset, timeInt, varargin)
%% DESCRIPTION
% 1. INTRODUCTION
% GETSPECTRA is used to retrieve scan from a dataset at a particular time
% or time interval. This function will update the finneeStc in case the
% extract spectra is to be save (optional) and the 2xm array spectraOut
% that contain the scan or sum of scans. This function work with profile or
% centroid spectrum, and extracted ions. 
%
% 2. PARAMETERS:
%   .required. GETSPECTRA requires at least 3 parameters
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
%       timeInt
%           Define the time interval to be used. It can be a single value
%           (closest time value) or a 2x1 array [tmin tmax]
%
%   .optionals. VARARGIN describes the optional paramters.  
%       'mzInt' follow by a 2x1 array
%           Allow to only take the data points with m/z are within the
%           defined interval
%       'noFig' 
%           No figures displayed
%       'indice'  
%           Work only with 'profile spectrum' and 'centroid spectrum'
%           dataset. If indice is used, timeInt is the scan number.
%       'reduced'
%           Remove null values in the spectra to decrease size
%
% 3. EXAMPLES:
%   spectraOut = getSpectra(finneeStc, 1, [10.1 10.5])
%
% 4. COPYRIGHT
% Copyright 2015-2016 G. Erny (guillaume@fe.up,pt), FEUP, Porto, Portugal

%% CORE OF THE FUNCTION
% 1. INITIALISATION
info.function.functionName =  'getSpectra';
info.function.description{1} = 'get the MS scan at a particular time or time interval';
info.function.matlabVersion = '8.5.0.197613 (R2015a)';
info.function.version = '14/01/2016';
info.function.ownerContact = 'guillaume@fe.up.pt';
[parameters, options] = ...
    initFunction(nargin, finneeStc, dataset, timeInt, varargin );
%INITFUNCTION used to verify the entries and load the optional and
% compulsory parameters



m = parameters.dataset;
fidReadDat = fopen(finneeStc.path2dat, 'rb');

%% FUNCTION CORE

axeX(:,1) = finneeStc.dataset{m}.axes.time.values;
if options.indice
    indTimeStt = parameters.xMin;
    indTimeEnd = parameters.xMax;
else
    indTimeStt = findCloser(parameters.xMin, axeX);
    indTimeEnd = findCloser(parameters.xMax, axeX);
end

spectraOut.title = ['MS spectra from ', num2str(axeX(indTimeStt)), ...
    ' to ', num2str(axeX(indTimeEnd)), ' ', ...
    finneeStc.dataset{m}.axes.time.unit];
spectraOut.data = [];

% 1. Checking the data type
switch finneeStc.dataset{m}.description.dataFormat
    case 'profile spectrum'
        spectraOut.plotType = 'profile';
        
        % 1. loading reference MZ axe
        index = finneeStc.dataset{m}.indexInDat(1, :);
        fseek(fidReadDat, index(1), 'bof');
        refMZ = ...
            fread(fidReadDat, [(index(2)-index(1))/(index(3)*4), index(3)], 'single');
        
        % 2. correcting MS spectra with inTimeStt
        index = finneeStc.dataset{m}.indexInDat(indTimeStt + 1, :);
        axeMZ = refMZ - index(4)*refMZ;
        ind2rem = axeMZ < parameters.mzMin | axeMZ > parameters.mzMax;
        axeMZ(ind2rem) = [];
        spectraOut.data(:,1)  = axeMZ;
        spectraOut.data(:,2) = 0;
        
        % 2. calculating the sum of MZ 
        for ii = indTimeStt:indTimeEnd
            index = finneeStc.dataset{m}.indexInDat(ii+1, :);
            fseek(fidReadDat, index(1), 'bof');
            
            switch index(6)
                case 2
                    MS = ...
                        fread(fidReadDat, [(index(2)-index(1))/(index(3)*2), index(3)], 'uint16');
                    
                case 4
                    MS = ...
                        fread(fidReadDat, [(index(2)-index(1))/(index(3)*4), index(3)], 'single');
                    
                case 8
                    MS = ...
                        fread(fidReadDat, [(index(2)-index(1))/(index(3)*8), index(3)], 'double');
            end
            
            MS(ind2rem) = [];
            
            spectraOut.data (:,2) = spectraOut.data(:,2) + MS;
                % NOTE: We assumed that the m/z axes in the same for all
                % scan this is not entirely true as there are some small
                % variation in the m/z values those seems small enough to
                % be discarded
                
        end
        
        if options.dtReduc
            % Data reduction for higher speed
            filterIn = true(length(spectraOut.data(:,1)), 1);
            filterIn(2:end-1) = spectraOut.data(1:end-2,2) == 0 & ...
                spectraOut.data(2:end-1,2) == 0 & ...
                spectraOut.data(3:end,2) == 0;
            spectraOut.data(filterIn, :) = [];
        end
        
    case 'centroid spectrum'
        spectraOut.plotType = 'stem';
        
        % 2.2 getting MS spectra within the interval indTimeStt:indTimeEnd
        for ii = indTimeStt:indTimeEnd
            index = finneeStc.dataset{m}.indexInDat(ii, :);
            fseek(fidReadDat, index(1), 'bof');
            
            switch index(6)
                case 2
                    MS = ...
                        fread(fidReadDat, [(index(2)-index(1))/(index(3)*2), index(3)], 'uint16');
                    
                case 4
                    MS = ...
                        fread(fidReadDat, [(index(2)-index(1))/(index(3)*4), index(3)], 'single');
                    
                case 8
                    MS = ...
                        fread(fidReadDat, [(index(2)-index(1))/(index(3)*8), index(3)], 'double');
            end
            
            ind2rem = MS(:,1) < parameters.mzMin |...
                MS(:,1) > parameters.mzMax;
            MS(ind2rem, :) = [];
            if isempty(spectraOut.data)
                spectraOut.data = MS;
            else
                [newSpectra, ~, ic] =...
                    unique([spectraOut.data(:,1); MS(:,1)]);
                newSpectra(:,2) = 0;
                fromMSSpectra = ic(1:length(spectraOut.data(:,1)));
                fromMS = ic(length(spectraOut.data(:,1))+1:end);
                newSpectra(fromMSSpectra,2) = ...
                    newSpectra(fromMSSpectra,2) + spectraOut.data(:,2);
                newSpectra(fromMS,2) = ...
                    newSpectra(fromMS,2) + MS(:,2);
                spectraOut.data = newSpectra;
            end
        end
    case 'ionic profile'
        spectraOut.plotType = 'stem';
        
        %2.3 getting each PIP
        for ii = 1:length(finneeStc.dataset{m}.indexInDat(:,1))
            index = finneeStc.dataset{m}.indexInDat(ii, :);
            fseek(fidReadDat, index(1), 'bof');
            
            switch index(6)
                case 2
                    PIP = ...
                        fread(fidReadDat, [(index(2)-index(1))/(index(3)*2), index(3)], 'uint16');
                    
                case 4
                    PIP = ...
                        fread(fidReadDat, [(index(2)-index(1))/(index(3)*4), index(3)], 'single');
                    
                case 8
                    PIP = ...
                        fread(fidReadDat, [(index(2)-index(1))/(index(3)*8), index(3)], 'double');
            end
            
            ind2rem = PIP(:,1) < indTimeStt | PIP(:,1) > indTimeEnd | ...
                PIP(:,3) < parameters.mzMin | PIP(:,3) > parameters.mzMax;
            PIP(ind2rem, :) = [];
            if ~isempty(PIP)
                MS = [PIP(:,3) PIP(:,2)];
                MS = sortrows(MS, 1);
                if length(MS(:,1)) ~= length(unique(MS(:,1)))
                    [newMS, ~, ic] = unique(MS(:,1));
                    newMS(:,2) = 0;
                    for jj = 1:length(ic)
                        newMS(ic(jj),2) = newMS(ic(jj),2) + MS(jj,2);
                    end
                    MS = newMS;
                end
                
                if isempty(spectraOut.data)
                    spectraOut.data = MS;
                else
                    [newSpectra, ~, ic] =...
                        unique([spectraOut.data(:,1); MS(:,1)]);
                    newSpectra(:,2) = 0;
                    fromMSSpectra = ic(1:length(spectraOut.data(:,1)));
                    fromMS = ic(length(spectraOut.data(:,1))+1:end);
                    newSpectra(fromMSSpectra,2) = ...
                        newSpectra(fromMSSpectra,2) + spectraOut.data(:,2);
                    newSpectra(fromMS,2) = ...
                        newSpectra(fromMS,2) + MS(:,2);
                    spectraOut.data = newSpectra;
                end
            end
        end
    otherwise
        DIL
end
fclose(fidReadDat);

spectraOut.axes.axeX.label = finneeStc.dataset{m}.axes.mz.label;
spectraOut.axes.axeX.unit = finneeStc.dataset{m}.axes.mz.unit;
spectraOut.axes.axeY.label = finneeStc.dataset{m}.axes.intensity.label;
spectraOut.axes.axeY.unit = finneeStc.dataset{m}.axes.intensity.unit;
spectraOut.info = info;

% 3. Plot if requested
if options.display
    switch spectraOut.plotType
        case 'profile'
            plot(spectraOut.data(:,1), spectraOut.data(:,2));
            title(spectraOut.title);
            xlabel([spectraOut.axes.axeX.label, ' / ', spectraOut.axes.axeX.unit]);
            ylabel([spectraOut.axes.axeY.label, ' / ', spectraOut.axes.axeY.unit]);
        case 'stem'
            stem(spectraOut.data(:,1), spectraOut.data(:,2), 'Marker', 'none');
            title(spectraOut.title);
            xlabel([spectraOut.axes.axeX.label, ' / ', spectraOut.axes.axeX.unit]);
            ylabel([spectraOut.axes.axeY.label, ' / ', spectraOut.axes.axeY.unit]);
        case 'bar'
            bar(spectraOut.data(:,1), spectraOut.data(:,2));
            title(spectraOut.title);
            xlabel([spectraOut.axes.axeX.label, ' / ', spectraOut.axes.axeX.unit]);
            ylabel([spectraOut.axes.axeY.label, ' / ', spectraOut.axes.axeY.unit]);
    end
end


%% NESTED FUNCTIONS
end

%% SUB FUNCTIONS
% 1. INITFUNCTION
% Function that get the input argument and check for errors
function [parameters, options] = ...
    initFunction(narginIn, finneeStc, dataset, timeInt, vararginIn )

options.display = 1;
options.indice = 0;
options.dtReduc = false;

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

if length(timeInt) == 1
    [parameters.xMin, parameters.xMax] = deal(timeInt);
elseif length(timeInt) == 2
    timeInt = sort(timeInt);
    parameters.xMin = timeInt(1);
    parameters.xMax = timeInt(2);
else
    error('line 253')
end

parameters.mzMin = 0;
parameters.mzMax = inf;

% 1.2. Check for options
if  narginIn > 3
    SFi = 1;
    while SFi <= length(vararginIn)
        switch vararginIn{SFi}
            case 'mzInt'
                if length(vararginIn{SFi+1}) == 2
                    mzInt = sort(vararginIn{SFi+1});
                    parameters.mzMin = mzInt(1);
                    parameters.mzMax = mzInt(2);
                end
                SFi = SFi + 2;
            case 'noFig'
                options.display = 0;
                SFi = SFi + 1;
            case 'reduced'
                options.dtReduc = true;
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
