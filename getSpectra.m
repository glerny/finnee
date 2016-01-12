function [spectraOut, finneeStc] = getSpectra(finneeStc, dataset, timeInt, varargin)
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
%           Allow to avoid displaying the resulting figure
%       'save2str'  
%           Will add the results in the list of trace associated with the 
%           dataset
%       'fid' followed by a file identifier
%           The file identifier shoud be related to the dat file, i.e.:
%           fid= fopen(finneeStc.dataset{m}.description.path2DatFile, 'br');
%           This function is to be used if GETSPECTRA is inside a loop to
%           avoid always openind and closing the dat file.
%
% 3. EXAMPLES:
%   spectraOut = getSpectra(finneeStc, 1, [10.1 10.5])
%
% 4. COPYRIGHT
% Copyright 2014-2015 G. Erny (guillaume@fe.up,pt), FEUP, Porto, Portugal

%% CORE OF THE FUNCTION
% 1. INITIALISATION
info.functionName = 'getSpectra';
info.description{1} = 'get the MS scan at a particular time or time interval';
info.matlabVersion = '8.5.0.197613 (R2015a)';
info.version = '25/06/2015_gle01';
info.ownerContact = 'guillaume@fe.up,pt';
spectraOut = [];
[parameters, options] = initFunction(nargin, finneeStc, dataset, timeInt, varargin );
%INITFUNCTION used to verify the entries and load the optional and
% compulsory parameters

m = parameters.dataset;
if options.openfile
    fidReadDat = fopen(finneeStc.dataset{m}.description.path2DatFile, 'rb');
else
    fidReadDat = options.fidReadDat;
    if fidReadDat <= 0
        fidReadDat = fopen(finneeStc.dataset{m}.description.path2DatFile, 'rb');
    end
end

%% FUNCTION CORE
index = finneeStc.dataset{m}.description.axe;
fseek(fidReadDat, index(1), 'bof');
axeX = fread(fidReadDat, [(index(2)-index(1))/(index(3)*8), index(3)], 'double');
if options.indice
    indTimeStt = parameters.xMin;
    indTimeEnd = parameters.xMax;
else
    indTimeStt = findCloser(parameters.xMin, axeX);
    indTimeEnd = findCloser(parameters.xMax, axeX);
end


% 1. Checking the data type
switch finneeStc.dataset{m}.description.dataFormat
    case 'profile spectrum'
        plotType = 'profile';
        
        %2.1 getting MS spectra within the interval indTimeStt:indTimeEnd
        for ii = indTimeStt:indTimeEnd
            index = finneeStc.dataset{m}.description.index2DotDat(ii, :);
            fseek(fidReadDat, index(1), 'bof');
            MS = ...
                fread(fidReadDat, [(index(2)-index(1))/(index(3)*8), index(3)], 'double');
            ind2rem = MS(:,1) < parameters.mzMin |...
                MS(:,1) > parameters.mzMax;
            MS(ind2rem, :) = [];
            
            if isempty(spectraOut)
                spectraOut = MS;
            else
                spectraOut(:,2) = spectraOut(:,2) + MS(:,2);
                % NOTE: We assumed that the m/z axes in the same for all
                % scan this is not entirely true as there are some small
                % variation in the m/z values those seems small enough to
                % be discarded
            end
        end
        
    case 'centroid spectrum'
        plotType = 'stem';        
        
        %2.2 getting MS spectra within the interval indTimeStt:indTimeEnd
        for ii = indTimeStt:indTimeEnd
            index = finneeStc.dataset{m}.description.index2DotDat(ii, :);
            fseek(fidReadDat, index(1), 'bof');
            MS = ...
                fread(fidReadDat, [(index(2)-index(1))/(index(3)*8), index(3)], 'double');
            % !Provisory solution to deal with different type of centroid
            % data
            typeCtr = length(MS(1,:));
            switch typeCtr
                case 2
                    ctrMaker = 'CompassXport';
                    ind2rem = MS(:,1) < parameters.mzMin |...
                        MS(:,1) > parameters.mzMax;
                    MS(ind2rem, :) = [];
                    if isempty(spectraOut)
                        spectraOut = MS;
                    else
                        [newSpectra, ~, ic] =...
                            unique([spectraOut(:,1); MS(:,1)]);
                        newSpectra(:,2) = 0;
                        fromMSSpectra = ic(1:length(spectraOut(:,1)));
                        fromMS = ic(length(spectraOut(:,1))+1:end);
                        newSpectra(fromMSSpectra,2) = ...
                            newSpectra(fromMSSpectra,2) + spectraOut(:,2);
                        newSpectra(fromMS,2) = ...
                            newSpectra(fromMS,2) + MS(:,2);
                        spectraOut = newSpectra;
                    end
                case 10
                    ctrMaker = 'gle1';
                    MS = [MS(:,7), MS(:,3)];
                    ind2rem = MS(:,1) < parameters.mzMin |...
                        MS(:,1) > parameters.mzMax;
                    MS(ind2rem, :) = [];
                    if isempty(spectraOut)
                        spectraOut = MS;
                    else
                        [newSpectra, ~, ic] =...
                            unique([spectraOut(:,1); MS(:,1)]);
                        newSpectra(:,2) = 0;
                        fromMSSpectra = ic(1:length(spectraOut(:,1)));
                        fromMS = ic(length(spectraOut(:,1))+1:end);
                        newSpectra(fromMSSpectra,2) = ...
                            newSpectra(fromMSSpectra,2) + spectraOut(:,2);
                        newSpectra(fromMS,2) = ...
                            newSpectra(fromMS,2) + MS(:,2);
                        spectraOut = newSpectra;
                    end
            end
        end
    case 'ionic profile'
        plotType = 'stem';
        
        %2.3 getting each PIP
        for ii = 1:length(finneeStc.dataset{m}.description.index2DotDat(:,1))
            index = finneeStc.dataset{m}.description.index2DotDat(ii, :);
            fseek(fidReadDat, index(1), 'bof');
            PIP = fread(fidReadDat, [(index(2)-index(1))/(index(3)*8), ...
                index(3)], 'double');
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
                
                if isempty(spectraOut)
                    spectraOut = MS;
                else
                    [newSpectra, ~, ic] =...
                        unique([spectraOut(:,1); MS(:,1)]);
                    newSpectra(:,2) = 0;
                    fromMSSpectra = ic(1:length(spectraOut(:,1)));
                    fromMS = ic(length(spectraOut(:,1))+1:end);
                    newSpectra(fromMSSpectra,2) = ...
                        newSpectra(fromMSSpectra,2) + spectraOut(:,2);
                    newSpectra(fromMS,2) = ...
                        newSpectra(fromMS,2) + MS(:,2);
                    spectraOut = newSpectra;
                end
            end
        end
    otherwise
        DIL
end
fclose(fidReadDat);

% 3. Plot and save if requested
if options.display
    if ~isempty(spectraOut)
    strName = ['MS spectra from ', num2str(parameters.xMin), ' to ', ...
        num2str(parameters.xMax), ' ', ...
        finneeStc.dataset{m}.description.timeUnit];
    switch plotType
        case 'profile'
            plot(spectraOut(:,1), spectraOut(:,2));
        case 'stem'
            stem(spectraOut(:,1), spectraOut(:,2), 'Marker', 'none');
        case 'barPlot'
            bar(spectraOut(:,1), spectraOut.traceOut(:,2), 1);  
    end
    title(strName);
    xlabel([finneeStc.dataset{m}.description.mzLabel,' / ',...
        finneeStc.dataset{m}.description.mzUnit]);
    ylabel([finneeStc.dataset{m}.description.intLabel, ' / ',...
        finneeStc.dataset{m}.description.intUnit]);
    end
end

if options.save2str
    finneeStc.dataset{m}.trace{end+1}.infoFunctionUsed.info = info;
    finneeStc.dataset{m}.trace{end}.infoFunctionUsed.parameters = parameters;
    finneeStc.dataset{m}.trace{end}.description.name = strName;
    finneeStc.dataset{m}.trace{end}.description.dateOfCreation = clock;
    finneeStc.dataset{m}.trace{end}.description.plotType = plotType;
    finneeStc.dataset{m}.trace{end}.description.axeX.label = ...
        finneeStc.dataset{m}.description.mzLabel;
    finneeStc.dataset{m}.trace{end}.description.axeX.unit = ...
        finneeStc.dataset{m}.description.mzUnit;
    finneeStc.dataset{m}.trace{end}.description.axeY.label = ...
        finneeStc.dataset{m}.description.intLabel;
    finneeStc.dataset{m}.trace{end}.description.axeY.unit = ...
        finneeStc.dataset{m}.description.intUnit;
    fidWriteTra = fopen( finneeStc.dataset{m}.description.path2DatFile, 'ab');
    fseek(fidWriteTra, 0,'eof');
    finneeStc.dataset{m}.trace{end}.index2DotDat  = ...
        [ftell(fidWriteTra), 0, 2];
    fwrite(fidWriteTra, [spectraOut(:,1) spectraOut(:,2)], 'double');
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
    initFunction(narginIn, finneeStc, dataset, timeInt, vararginIn )

options.display = 1;
options.save2str = 0;
options.openfile = 1;
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
            case 'save2str'
                options.save2str = 1;
                SFi = SFi + 1;
            case 'fid'
                options.openfile = 0;
                options.fidReadDat = vararginIn{SFi+1};
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
