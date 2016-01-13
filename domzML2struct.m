function finneeStc = domzML2struct(varargin)
%% DESCRIPTION 
% 1. Introduction
% DOMZML2STRUCT transforms a mzML file to a a matlab structure and an associated data file. This function has been
% specifically design to handle datasets obtained using separation techniques
% hyphenated with high resolution mass spectrometers. Before using DOMZML2STRUCT the raw datafile need to be converted to mzML [1]. The mzML format is a
% markup text language supported by the human proteome organization [2]
%  that promote a standardized format for mass spectrometry data. Conversion of a
% raw dataset to mzML file can be made via freeware such as msConvert from
% ProteoWizard [3] or by tools made available by the instruments
% manufacturers [4]. DOMZML2STRUCT is suitable for profiles and centroid spectrum datasets.
%
% 2. Usage
% *finneeStc* is a Matlab(R) structure used to record all
% information link to a specific dataset. This structure has been design in
% order to store original information and to allow additional
% transformations *finneeStc* contain multiple fields organised in
% substructures. The initial level will contain at least three main
% substructures:
%
% * *_finneeStc.infoFunctionUsed_* This is the substructure used to record
% information about the function used to convert the mzML files to the
% Matlab(R) structure
% * *_finneeStc.infoRun_* This is the substructure used to record all the
% information contained in the mzML file up to the mzML tag '<run>'. This substructure copy the mzML format to allow for easy conversion.
% * *_finneeStc.dataset{n}_* are the substructures were  the various scans are
% indexed. a finneeStc can contain multiple datasets, for exemple a
% dataset containing the profile spectra, one containing the centroid
% spectra and one containing the ionic profiles. Those substructures are
% divided in three substructures.
%
% * *_finneeStc.dataset{n}.infoFunctionUsed_* contain the information about
% the function used to obtain this dataset.
% * *_finneeStc.dataset{n}.description_*. This is the key structure to
% retrieve each scans. Those are store in a secondary binary files with the
% extension .dat. Fields in desciption are:
%%
%  path2DatFile_  : full address to the binary data file
%  mzStart_       :
%  mzEnd_         :
%  intMin_        : 
%  intMax_        : 
%  timeStart_     : 
%  timeEnd_       :
%  dataFormat_    : 'centroid spectrum', 'profile spectrum' or 'ionic
% profile'
%  timeLabel_     : i.e. Time
%  timeUnit_      : min
%  mzLabel_       : 
%  mzUnit_        :
%  intLabel_      :
%  intUnit_       :
%  index2DotDat_  : this is a yx3 array where y is the number of scans
% recorded during the run. This allow to obtain rapidly any recorded scans.
% for example to obtain the 50th scans the code will be
%%
% 
%   fileDat = finneeStc.dataset{1}.description.mzStart
%   fidReadDat = fopen(fileDat, 'rb');
%   index = finneeStc.dataset{1}.description.index2DotDat(50, :)
%   fseek(fidReadDat, ind(1), 'bof');
%   MSScan = fread(fidReadDat, [(ind(2)-ind(1))/(8*ind(1)) ind(3)], 'double');
%   fclose(fidReadDat)
%
%%
%  axe_ : 1x3 array with the sane rule as _index2DotDat_, allows to
% recover the time axe.
%%
% * *_finneeStc.trace{m}_* are structures where every traces (i.e. any 2
% dimensionals represenation such as chromatogram, MS spectra at a given
% time or intervals,...) are recorded.
%
%% 
% *3. Optional parameters.*
%%
% varagin is used to indicate optionals parameters. The following
% parameters are recognized  
%%
%  'tMin:Max'  followed by a 2x1 array of integer to set the  time interval (default is (0 inf))
%  'MZMin:Max' followed by a 2x1 array to set the MZ interval (default is (0 inf))
%  'display'   followed by 'Off' remove all display text in the Matlab command window.
%  'overwrite' Overwrite potentially existing files
%
%%  EXAMPLES:
%   finneeStc = mzML2struct('tMin:Max', [5 25], 'display', 'Off', 'overwrite');
%
%% REFERENCES
% [1] 
% [2] <http://www.psidev.info/mzml_1_0_0%20>
% [3] <http://proteowizard.sourceforge.net/>
% [4] <http://tools.proteomecenter.org/wiki/index.php?title=Formats:mzXML>

%% ACKNOWLEDGEMENTS
%% COPYRIGHT
% Copyright 2014-2016 G. Erny (guillaume@fe.up,pt), FEUP, Porto, Portugal

%% CORE OF THE FUNCTION
% 1. INITIALISATION
info.functionName = 'MSdata2struct';
info.description{1} = 'Main function';
info.description{2} = 'Used to convert mzML datafiles to Matlab structures';
info.matlabVersion = '8.5.0.197613 (R2015a)';
info.version = '25/06/2015_gle01';
info.ownerContact = 'guillaume.erny@finnee.com';
finneeStc.infoFunctionUsed.info = info;
[finneeStc.infoFunctionUsed.parameters, options] = ...
    initFunction(nargin, varargin);
finneeStc.infoFunctionUsed.errors = {};
% INITFUNCTION is used to test the entries, load the target MS dataset and
% load the default values

originalFile = finneeStc.infoFunctionUsed.parameters.fileIn;
fileID = finneeStc.infoFunctionUsed.parameters.fileID;
folderOut = finneeStc.infoFunctionUsed.parameters.folderOut;
try
    fidRead = fopen(originalFile, 'r'); % orifinal mzML file
    fidWriteDat = fopen(fullfile(folderOut, [fileID, '.dat']), 'wb'); % new file that will receive the binary data
catch
    error('myApp:argChk', ...
        'Error while opening the files. Type help hyphMSdata2struct for more information')
end
if fidRead == -1 && fidWriteDat == -1 
    error('myApp:argChk', ...
        'Error while opening the files. Type help hyphMSdata2struct for more information')
end

% 2. RECORD GENERAL INFORMATIONS ABOUT THE RUN
frewind(fidRead)
curLine = fgetl(fidRead);
[LDR, curLine] = getMZMLCamp(curLine, fidRead);
while ~strcmp(LDR.label, 'run')
    if strcmp(LDR.label, '?xml')
        finneeStc.infoRun.xml.attributes = LDR.attributes;
        finneeStc.infoRun.xml.text = LDR.text;
    elseif strcmp(LDR.label, 'mzML')
        finneeStc.infoRun.mzML.attributes = LDR.attributes;
        finneeStc.infoRun.mzML.text = LDR.text;
    else
        s = addChildNode(LDR);
        finneeStc.infoRun.(LDR.label) = s.(LDR.label);
    end
    curLine = fgetl(fidRead);
    if ~ischar(curLine)
         error('myApp:endOfFile', ...
             'The mzML file is not complete')
    end
    [LDR, curLine] = getMZMLCamp(curLine, fidRead);
end
finneeStc.infoRun.run.attributes = LDR.attributes;
finneeStc.infoRun.run.text = LDR.text;

% 3. CHECK FIRTH SCAN FOR MISSING INOFRMATION
finneeStc.dataset{1}.infoFunctionUsed.info = info;
finneeStc.dataset{1}.infoFunctionUsed.parameters = ...
    finneeStc.infoFunctionUsed.parameters;
errors = {};
finneeStc.dataset{1}.infoFunctionUsed.errors = {};

% Get dataset axesLabel and axesUnit by processing the first scan
% Initialisation and default values
timeLabel = 'Time'; timeUnit = ''; mzLabel = 'Mass';
mzUnit = ''; intLabel = 'Intensity'; intUnit = '';
[mzMin, intMin] = deal(inf); [mzMax, intMax] = deal(0);  
axeX = []; TICP = [];  BPP = []; mzBPP = []; MSIndex = []; data2keep = [];
compression = ''; dataFormat = ''; 
decimals = finneeStc.dataset{1}.infoFunctionUsed.parameters.decimals;
datasetType = finneeStc.infoRun.fileDescription.fileContent{1}.cvParam{2}.attributes{3}.field{1};
curLine = fgetl(fidRead);
[LDR, curLine] = getMZMLCamp(curLine, fidRead);
if options.display
    disp('processing scan 1 out of ???')
end
boolArray = 1;
while ~strcmp(LDR.label, '/spectrum')
    if strcmp(LDR.label, 'spectrumList')
        [bool, ~, provField] = fieldfind( LDR.attributes, 'count');
        if bool
            scanCount = str2double(provField{1});
        else
            scanCount = 0;
            errors{end+1} = 'attribute count not present in tag spectrumList';
        end
    elseif strcmp(LDR.label, 'cvParam')
        [bool, ~, field] = fieldfind( LDR.attributes, 'accession');
        if bool
            switch field{1}
                case 'MS:1000016'
                    [~, ~, provField] = fieldfind( LDR.attributes, 'value');
                    axeX(end+1) = str2double(provField{1});
                    [~, ~, provField] = fieldfind( LDR.attributes, 'unitName');
                    timeUnit = provField{1};
                case 'MS:1000574'
                    compression = field{1};
                case 'MS:1000514'
                    [~, ~, provField] = fieldfind( LDR.attributes, 'unitName');
                    mzUnit = provField{1};
                case 'MS:1000523'
                    dataFormat = field{1};
                case 'MS:1000515'
                    [~, ~, provField] = fieldfind( LDR.attributes, 'unitName');
                    intUnit = provField{1};
            end
        end
    elseif strcmp(LDR.label, 'binary')
        if axeX(end) >= finneeStc.dataset{1}.infoFunctionUsed.parameters.xMin && ...
                axeX(end) <= finneeStc.dataset{1}.infoFunctionUsed.parameters.xMax
            data2keep(end+1) = 1;
            input = LDR.text;
            switch dataFormat
                case 'MS:1000523'
                    output = base64decode(input);
                otherwise
                    error('precision not recognized')
            end
            switch compression
                case 'MS:1000574'
                    output = zlibdecode(output);
                otherwise
                    error('compression not recognized')
            end
            output= typecast(uint8(output),'double');
            if boolArray
                boolArray = 0;
                mzValue = output;
            else
                boolArray = 1;
                intValue = output;
                ind2rem = mzValue < finneeStc.dataset{1}.infoFunctionUsed.parameters.mzMin | ...
                    mzValue > finneeStc.dataset{1}.infoFunctionUsed.parameters.mzMax;
                mzValue(ind2rem) = [];
                intValue(ind2rem) = [];
                if isnan(decimals) && strcmp(datasetType, 'centroid spectrum')
                    answer = inputdlg('Enter the number of decimals to keep in the accurate mass', ...
                        'decimals', 1, {'5'});
                    if isempty(answer)
                        error('myApp:argChk', 'User cancel')
                    elseif isempty(answer{1})
                        error('myApp:argChk', 'User cancel')
                    elseif str2double(answer{1}) < 1
                        error('myApp:argChk', 'invalid value for the number of decimals. Should be bigger or equal than 1')

                    end
                    decimals = round(str2double(answer{1}));
                    finneeStc.dataset{1}.infoFunctionUsed.parameters.decimals = decimals;
                end
                if strcmp(datasetType, 'centroid spectrum')
                    mzValue = round(mzValue*10^decimals)/10^decimals;
                end
                
                % calculates profiles and limits
                if mzMin > min(mzValue), mzMin = min(mzValue); end
                if mzMax < max(mzValue), mzMax = max(mzValue); end
                if intMin > min(intValue), intMin = min(intValue); end
                if intMax < max(intValue), intMax = max(intValue); end
                TICP(end+1) = sum(intValue);
                [BPP(end+1), indMax] = max(intValue);
                mzBPP(end+1) = mzValue(indMax);
                
                % write 
                posIni =  ftell(fidWriteDat);
                fwrite(fidWriteDat, ...
                    [mzValue intValue], 'double');
                MSIndex = [MSIndex; [posIni, ftell(fidWriteDat), 2]];
            end
        else
            data2keep(end+1) = 0;
        end
                        
    end
    curLine = fgetl(fidRead);
    if ~ischar(curLine)
         error('myApp:endOfFile', ...
             'The mzML file is not complete')
    end
    [LDR, curLine] = getMZMLCamp(curLine, fidRead);
end

% 4. DO FOR EACH SCAN
while ~strcmp(LDR.label, '/run')
    if strcmp(LDR.label, 'spectrum'), boolArray = 1; end
    if strcmp(LDR.label, 'cvParam') && ~isempty(LDR.attributes)
        [~, ~, field] = fieldfind( LDR.attributes, 'accession');
        if strcmp(field{1}, 'MS:1000016')
            [~, ~, provField] = fieldfind( LDR.attributes, 'value');
            axeX(end+1) = str2double(provField{1});
            if options.display
                disp(['processing scan ', num2str(length(axeX)), ...
                    ' out of ', num2str(scanCount)])
            end
        end
    end
    if strcmp(LDR.label, 'binary')
        if axeX(end) >= finneeStc.dataset{1}.infoFunctionUsed.parameters.xMin && ...
                axeX(end) <= finneeStc.dataset{1}.infoFunctionUsed.parameters.xMax
            data2keep(end+1) = 1;
            input = LDR.text;
            switch dataFormat
                case 'MS:1000523'
                    output = base64decode(input);
                otherwise
                    error('precision not recognized')
            end
            switch compression
                case 'MS:1000574'
                    output = zlibdecode(output);
                otherwise
                    error('compression not recognized')
            end
            output= typecast(uint8(output),'double');
            if boolArray
                boolArray = 0;
                mzValue = output;
            else
                boolArray = 1;
                intValue = output;
                ind2rem = mzValue < finneeStc.dataset{1}.infoFunctionUsed.parameters.mzMin | ...
                    mzValue > finneeStc.dataset{1}.infoFunctionUsed.parameters.mzMax;
                mzValue(ind2rem) = [];
                intValue(ind2rem) = [];
                if strcmp(datasetType, 'centroid spectrum')
                    mzValue = round(mzValue*10^decimals)/10^decimals;
                end
                
                % calculates profiles and limits
                if mzMin > min(mzValue), mzMin = min(mzValue); end
                if mzMax < max(mzValue), mzMax = max(mzValue); end
                if intMin > min(intValue), intMin = min(intValue); end
                if intMax < max(intValue), intMax = max(intValue); end
                TICP(end+1) = sum(intValue);
                [BPP(end+1), indMax] = max(intValue);
                mzBPP(end+1) = mzValue(indMax);
                % write
                posIni =  ftell(fidWriteDat);
                fwrite(fidWriteDat, ...
                    [mzValue intValue], 'double');
                MSIndex = [MSIndex; [posIni, ftell(fidWriteDat), 2]];
            end
        else
            data2keep(end+1) = 0;
        end
    end
    
    curLine = fgetl(fidRead);
    if ~ischar(curLine)
        error('myApp:endOfFile', ...
            'The mzML file is not complete')
    end
    [LDR, curLine] = getMZMLCamp(curLine, fidRead);
end
axeX(data2keep == 0) = [];
fclose(fidRead);

% 5. CONCLUSIONS
finneeStc.dataset{1}.infoFunctionUsed.errors = errors;
save2struc()
save(fullfile(finneeStc.infoFunctionUsed.parameters.folderOut, ...
    [finneeStc.infoFunctionUsed.parameters.fileID '.fin']), 'finneeStc', '-mat')
fclose(fidWriteDat);

%% SUB FUNCTIONS AND NESTED FUNCTIONS
% 1. NESTED-FUNCTION: outStrc = addChildNode(LDR)
    function outStrc = addChildNode(LDR)
        outStrc = struct;
        outStrc.(LDR.label).attributes = LDR.attributes;
        outStrc.(LDR.label).text= LDR.text;
        
        if LDR.open
            endLabel = ['/', LDR.label];
            parentLabel = LDR.label;
            while 1
                curLine = fgetl(fidRead);
                [LDR, ~] = getMZMLCamp(curLine, fidRead);
                LDR.label
                if strcmp(endLabel, LDR.label)
                    break
                end
                if ~isfield(outStrc.(parentLabel), LDR.label)
                    outStrc.(parentLabel).(LDR.label) = {};
                end
                s = addChildNode(LDR);
                outStrc.(parentLabel).(LDR.label){end+1} = s.(LDR.label);
            end
        end
        
    end

% 2. NESTED-FUNCTION:  save2struc()
    function save2struc()
    % Save description of dataset
    finneeStc.dataset{1}.description.path2DatFile = fullfile(folderOut, [fileID, '.dat']);
    finneeStc.dataset{1}.description.mzStart = mzMin;
    finneeStc.dataset{1}.description.mzEnd = mzMax;
    finneeStc.dataset{1}.description.intMin = intMin;
    finneeStc.dataset{1}.description.intMax = intMax;
    finneeStc.dataset{1}.description.timeStart = axeX(1);
    finneeStc.dataset{1}.description.timeEnd = axeX(length(axeX));
    finneeStc.dataset{1}.description.dataFormat = datasetType;
    finneeStc.dataset{1}.description.index2DotDat = MSIndex;
    
    % Record axeX in the *dat* file the rest in the *tra* file
    finneeStc.dataset{1}.description.axe(1) = ftell(fidWriteDat);
    fwrite(fidWriteDat, axeX, 'double');
    finneeStc.dataset{1}.description.axe(2) = ftell(fidWriteDat);
    finneeStc.dataset{1}.description.axe(3) = 1;
    finneeStc.dataset{1}.description.timeLabel = timeLabel;
    finneeStc.dataset{1}.description.timeUnit = timeUnit;
    finneeStc.dataset{1}.description.mzLabel = mzLabel;
    finneeStc.dataset{1}.description.mzUnit = mzUnit;
    finneeStc.dataset{1}.description.intLabel = intLabel;
    finneeStc.dataset{1}.description.intUnit = intUnit;
    
    % Record profiles
    % ** TICP
    finneeStc.dataset{1}.trace{1}.infoFunctionUsed.info = info;
    finneeStc.dataset{1}.trace{1}.infoFunctionUsed.parameters = {};
    finneeStc.dataset{1}.trace{1}.description.name = 'Total Ion Current Profile (dataset 1)';
    finneeStc.dataset{1}.trace{1}.description.dateOfCreation = clock;
    finneeStc.dataset{1}.trace{1}.description.plotType = 'profile';
    finneeStc.dataset{1}.trace{1}.description.axeX.label = timeLabel;
    finneeStc.dataset{1}.trace{1}.description.axeX.unit = timeUnit;
    finneeStc.dataset{1}.trace{1}.description.axeY.label = intLabel;
    finneeStc.dataset{1}.trace{1}.description.axeY.unit = intUnit;
    finneeStc.dataset{1}.trace{1}.index2DotDat  = [ftell(fidWriteDat), 0, 2];
    fwrite(fidWriteDat, [axeX TICP], 'double');
    finneeStc.dataset{1}.trace{1}.index2DotDat(2) = ftell(fidWriteDat);
    
    % ** BPP
    finneeStc.dataset{1}.trace{2}.infoFunctionUsed.info = info;
    finneeStc.dataset{1}.trace{2}.infoFunctionUsed.parameters = {};
    finneeStc.dataset{1}.trace{2}.description.name = 'Base Peak Profile (dataset 1)';
    finneeStc.dataset{1}.trace{2}.description.dateOfCreation = clock;
    finneeStc.dataset{1}.trace{2}.description.plotType = 'profile';
    finneeStc.dataset{1}.trace{2}.description.axeX.label = timeLabel;
    finneeStc.dataset{1}.trace{2}.description.axeX.unit = timeUnit;
    finneeStc.dataset{1}.trace{2}.description.axeY.label = intLabel;
    finneeStc.dataset{1}.trace{2}.description.axeY.unit = intUnit;
    finneeStc.dataset{1}.trace{2}.index2DotDat  = [ftell(fidWriteDat), 0, 2];
    fwrite(fidWriteDat, [axeX BPP], 'double');
    finneeStc.dataset{1}.trace{2}.index2DotDat(2) = ftell(fidWriteDat);
    
    % ** mzBPP
    finneeStc.dataset{1}.trace{3}.infoFunctionUsed.info = info;
    finneeStc.dataset{1}.trace{3}.infoFunctionUsed.parameters = {};
    finneeStc.dataset{1}.trace{3}.description.name = 'm/z @ Base Peak (dataset 1)';
    finneeStc.dataset{1}.trace{3}.description.dateOfCreation = clock;
    finneeStc.dataset{1}.trace{3}.description.plotType = 'profile';
    finneeStc.dataset{1}.trace{3}.description.axeX.label = timeLabel;
    finneeStc.dataset{1}.trace{3}.description.axeX.unit = timeUnit;
    finneeStc.dataset{1}.trace{3}.description.axeY.label = mzLabel;
    finneeStc.dataset{1}.trace{3}.description.axeY.unit = mzUnit;
    finneeStc.dataset{1}.trace{3}.index2DotDat  = [ftell(fidWriteDat), 0, 2];
    fwrite(fidWriteDat, [axeX mzBPP], 'double');
    finneeStc.dataset{1}.trace{3}.index2DotDat(2) = ftell(fidWriteDat);
    end
end

% 3. SUB-FUNCTION: [parameters, options] = initFunction(narginIn, vararginIn)
function [parameters, options] = ...
    initFunction(narginIn, vararginIn)

% Check for obligatory parameters
if narginIn < 0 % check the number of input parameters
    error('myApp:argChk', ...
        ['Wrong number of input arguments. \n', ...
        'Type help MSdata2struct for more information']);
end

% Load default parameters and options
ext = 'pwd/*.mzML';
txtStg = 'Select the mzML file to load';
parameters.xMin = 0;
parameters.xMax = inf;
parameters.mzMin = 0;
parameters.mzMax = inf;
parameters.decimals = nan;
options.display = 1;
options.ifExist = 'Stop';

% Decipher varargin
if  narginIn > 0
    SFi = 1;
    length(vararginIn)
    while SFi <= length(vararginIn)
        switch vararginIn{SFi}
            case 'tMin:Max'
                if isnumeric(vararginIn{SFi+1}) 
                    if  vararginIn{SFi+1}(2) >  vararginIn{SFi+1}(1) > 0;
                        parameters.xMin = vararginIn{SFi+1}(1);
                        parameters.xMax = vararginIn{SFi+1}(2);
                    else
                         error('myApp:argChk', ...
                             ['incorrect set values for ''tMin:Max''',...
                             ' the condition tMax>tMin>0 is not verified'])
                    end
                else
                    error('myApp:argChk', ...
                             ['incorrect set values for ''tMin:Max''',...
                             ' the following parameter should be a 2x1 array of numbers'])
                end 
                SFi = SFi +2;
            case 'MZMin:Max'
                if isnumeric(vararginIn{SFi+1}) 
                    if  vararginIn{SFi+1}(2) >  vararginIn{SFi+1}(1) > 0;
                        parameters.mzMin = vararginIn{SFi+1}(1);
                        parameters.mzMax = vararginIn{SFi+1}(2);
                        SFi = SFi +2;
                    else
                         error('myApp:argChk', ...
                             ['incorrect set values for ''MZMin:Max''',...
                             ' the condition MZMax>MZMin>0 is not verified'])
                    end
                else
                    error('myApp:argChk', ...
                             ['incorrect set values for ''MZMin:Max''',...
                             ' the following parameter should be a 2x1 array of numbers'])
                end 
            case 'display'
                if strcmpi(vararginIn{SFi+1}, 'off')
                    options.display = 0;
                elseif strcmpi(vararginIn{SFi+1}, 'on')
                    options.display = 1;
                else
                    error('myApp:argChk', ...
                        ['incorrect argument for ''display''',...
                        ' the following parameter should either be ''on'' or ''off'''])
                end
                SFi = SFi +2;
            case 'overwrite'
                options.ifExist = 'Overwrite';
                SFi = SFi + 1;
            case 'decimals'
                if isnumeric(vararginIn{SFi+1})
                    if vararginIn{SFi+1} >= 1
                        parameters.decimals = int32(vararginIn{SFi+1});
                    else
                        error('myApp:argChk', ...
                            ['incorrect argument for ''decimals''',...
                            ' the following paramters should be >= 1'])
                    end
                else
                      error('myApp:argChk', ...
                            ['incorrect argument for ''decimals''',...
                            ' the following paramters should be numeric'])
                end
                SFi = SFi +1;
            otherwise
                error('myApp:argChk', ...
                    [vararginIn{SFi} ' is not a recognized PropertyName'])
        end
    end
end

% Ask and load files
[fileName, pathName] = uigetfile(ext, txtStg);
if ~ischar(fileName) && ~ischar(pathName)
    error('myApp:argChk', 'User cancel file selection');
end
parameters.fileIn = fullfile(pathName, fileName);
folderOut = uigetdir(pwd, 'Select the folder of destination');
if ~ischar(folderOut)
    error('myApp:argChk', 'User cancel');
end
parameters.folderOut = folderOut;
[~, dfltName, ~] = fileparts(parameters.fileIn);
answer = inputdlg('Enter a generic name for the file', ...
    'fileID', 1, {dfltName});
if isempty(answer)
    error('myApp:argChk', 'User cancel')
elseif isempty(answer{1})
    error('myApp:argChk', 'User cancel')
end
parameters.fileID = answer{1};

% If target files already exist
if strcmp(options.ifExist, 'Stop')
    if exist(fullfile(folderOut, [parameters.fileID '.fin']), 'file') == 2 || ...
            exist(fullfile(folderOut, [parameters.fileID '.tra']), 'file') == 2 || ...
            exist(fullfile(folderOut, [parameters.fileID '.dat']), 'file') == 2
        error('myApp:argChk', ...
            'file(s) of destination already exist, change name or select overwrite')
    end
end
end

% 4. SUB-FUNCTION [LDR, curLine] = getMZMLCamp(curLine, fidRead)
function [LDR, curLine] = getMZMLCamp(curLine, fidRead)
% Get the data_label and associated data-set from mzML line(s)
% LDR is a structure with 4 files: Labels (xml tag), attributes (a
% structure with multiple pair label/filed), Text and open (1 if the filed
% have child nodes.
curLine = strtrim(curLine); %remove leading and trailing white spece
if ~strcmp(curLine(1), '<')
    error('myApp:argChk', 'Incorrect data')
end

% find all opening bracket
LDR.label = '';
LDR.attributes = {};
LDR.text = {};
LDR.open = 1;

% find all opening bracket
indSt = strfind(curLine, '<');
if length(indSt) == 1
    
    % No intercalling text
    if strcmp(curLine(end), '>') %classical case (finish)
        if curLine(end-1) == '/'
             str2decode = curLine(2:end-2);
             LDR.open = 0;
        else
             str2decode = curLine(2:end-1);
             LDR.open = 1;
        end
        ind2att = strfind(str2decode, '="');
        if isempty(ind2att)
            LDR.label = str2decode;
        else
            ind2blk = strfind(str2decode, ' ');
            ind2dbc =  strfind(str2decode, '"');
            LDR.label = str2decode(1:ind2blk(1)-1);
            for ii = 1:length(ind2att)
                istt = find(ind2blk < ind2att(ii), 1, 'last');
                iend = find(ind2dbc > ind2att(ii), 1, 'first');
                istt = ind2blk(istt);
                iend = ind2dbc(iend+1);
                cutstr = str2decode(istt+1:iend-1);
                indEqu = strfind(cutstr, '=');
                LDR.attributes{end+1}.label = cutstr(1:indEqu(1)-1);
                LDR.attributes{end}.field{1} = cutstr(indEqu(1)+2:end);
            end
        end
        
    else        %case where field as entry in it
        str2decode = curLine(2:end);
        ind2att = strfind(str2decode, '="');
        if isempty(ind2att)
            LDR.label = str2decode;
        else
            ind2blk = strfind(str2decode, ' ');
            LDR.label = str2decode(1:ind2blk(1)-1);
            for ii = 1:length(ind2att)-1
                istt = find(ind2blk < ind2att(ii), 1, 'last');
                iend = find(ind2blk > ind2att(ii), 1, 'first');
                istt = ind2blk(istt);
                if isempty(iend)
                    iend = length(cutstr);
                else
                    iend = ind2blk(iend);
                end
                cutstr = str2decode(istt+1:iend-1);
                indEqu = strfind(cutstr, '=');
                LDR.attributes{end+1}.label = cutstr(1:indEqu(1)-1);
                indGlm = strfind(cutstr, '"');
                LDR.attributes{end}.field{1} = ...
                    cutstr(indGlm(1)+1:indGlm(2)-1);
            end
            ii = ii +1;
            istt = find(ind2blk < ind2att(ii), 1, 'last');
            istt = ind2blk(istt);
            iend = length(str2decode);
            cutstr = str2decode(istt+1:iend-1);
            indEqu = strfind(cutstr, '=');
            LDR.attributes{end+1}.label = cutstr(1:indEqu(1)-1);
            indGlm = strfind(cutstr, '"');
            LDR.attributes{end}.field{1} = ...
                cutstr(indGlm(1)+1:end);
            while 1
                curLine = fgetl(fidRead);
                if strcmp(curLine(end-1:end), '/>')
                    LDR.attributes{end}.field{end+1} = curLine(1:end-2);
                    LDR.open = 0;
                    break
                else
                    LDR.attributes{end}.field{end+1} = curLine;
                    LDR.open = 1;
                end
            end
        end
        
    end
elseif length(indSt) == 2 % More than two starting bracket, check for text
    indEn = strfind(curLine, '>');
    if strcmp(curLine(indSt(1)+1:indEn(1)-1), 'binary')
        LDR.label = 'binary';
        LDR.attributes = {};
        LDR.text = curLine(indEn(1)+1:indSt(2)-1);
    else
        SHOULDNOTHAPPENWITHmzMLFILES
    end
    
else
    SHOULDNOTHAPPENWITHmzMLFILES
end
end
