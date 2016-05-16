function finneeStc = domzML2struct(varargin)
%% DESCRIPTION 
% For detail information and tutorials, go to
% https://finneeblog.wordpress.com/. This work is still in progress and is
% not bug free. For questions, help, suggestions or proposition of
% collaboration please contact me via finneeblog or by email
% (guillaume@fe.up,pt).
%
% 1. INTRODUCTION 
% DOMZML2STRUCT transforms a mzML file to a a matlab
% structure and an associated data file. This function has been design to
% handle large datasets obtained using separation techniques hyphenated
% with high resolution mass spectrometers. Before using DOMZML2STRUCT the
% original proprietary file need to be converted to the mzML format.
% DOMZML2STRUCT is suitable for profiles and centroid spectrum datasets.
%
% 2. INPUT PARAMETERS
%   .required
%       None, the input file, destination files and folder are selected
%       within the function
%
%   .optional VARARGIN acepts the following values
%       'dispOff'   remove all display text in the Matlab command window.
%       'overwrite' Overwrite potentially existing files
%
% 3. OUTPUT PARAMETERS
%   .finneeStc 
%       is the structure that contain all information that was contain in
%       the mzML file.
%
% 4.  EXAMPLES
%   finneeStc = mzML2struct('tMin:Max', [5 25], 'dispOff', 'overwrite');
%
% 5. COPYRIGHT
% Copyright 2015-2016 G. Erny (guillaume@fe.up,pt), FEUP, Porto, Portugal

%% CORE OF THE FUNCTION
% 1. INITIALISATION
info.function.functionName = 'domzML2struct';
info.function.description{1} = 'Main function';
info.function.description{2} = 'Used to convert mzML datafiles to Matlab structures';
info.function.matlabVersion = '8.5.0.197613 (R2015a)';
info.function.version = '12/01/2016';
info.function.ownerContact = 'guillaume.erny@finnee.com';
finneeStc.info = info;
warning('off','all')

% INITFUNCTION is used to test the entries, load the target MS dataset and
% load the default values
[finneeStc.info.parameters, options] = initFunction(nargin, varargin);

finneeStc.info.errors = {};
originalFile = finneeStc.info.parameters.fileIn;
fileID = finneeStc.info.parameters.fileID;
folderOut = finneeStc.info.parameters.folderOut;
finneeStc.path2dat = fullfile(folderOut, [fileID, '.dat']);
try
    fidRead = fopen(originalFile, 'r'); % original mzML file
    fidWriteDat = fopen(finneeStc.path2dat, 'wb'); % new file that will store data
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
        finneeStc.mzML.xml.attributes = LDR.attributes;
        finneeStc.mzML.xml.text = LDR.text;
    elseif strcmp(LDR.label, 'indexedmzML')
        finneeStc.mzML.indexedmzML.attributes = LDR.attributes;
        finneeStc.mzML.indexedmzML.text = LDR.text;
    elseif strcmp(LDR.label, 'mzML')
        finneeStc.mzML.mzML.attributes = LDR.attributes;
        finneeStc.mzML.mzML.text = LDR.text;
    else
        s = addChildNode(LDR);
        finneeStc.mzML.(LDR.label) = s.(LDR.label);
    end
    curLine = fgetl(fidRead);
    if ~ischar(curLine)
         error('myApp:endOfFile', ...
             'The mzML file is not complete')
    end
    [LDR, curLine] = getMZMLCamp(curLine, fidRead);
    assignin('base', 'mzML', finneeStc.mzML)
end
finneeStc.mzML.run.attributes = LDR.attributes;
finneeStc.mzML.run.text = LDR.text;

% 3. CHECK FIRST SCAN FOR MISSING INOFRMATION

% Get dataset axesLabel and axesUnit by processing the first scan
% Initialisation and default values
timeLabel = 'Time'; timeUnit = ''; mzLabel = 'Mass';
mzUnit = ''; intLabel = 'Intensity'; intUnit = '';
[mzMin, intMin] = deal(inf); [mzMax, intMax] = deal(0);  
axeX = []; TICP = [];  BPP = []; mzBPP = []; MSIndex = []; data2keep = [];
compression = ''; dataFormat = ''; 
keptPl = ftell(fidRead)
curLine = fgetl(fidRead)
[LDR, curLine] = getMZMLCamp(curLine, fidRead);
if options.display
    disp('processing scan 1 out of ???')
end
boolArray = 1;

% MODIDICATION 16/05/2016 FOR  datasetType = 'profile spectrum';
% SCAN WHOLE DATASET TO GENERATE mzAxe 
%record position where run start

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
                case 'MS:1000127'
                    datasetType = 'centroid spectrum';
                case 'MS:1000128'
                    datasetType = 'profile spectrum';
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
    end
    curLine = fgetl(fidRead);
    if ~ischar(curLine)
         error('myApp:endOfFile', ...
             'The mzML file is not complete')
    end
    [LDR, curLine] = getMZMLCamp(curLine, fidRead);
end

fseek(fidRead, keptPl, 'bof') % Go back at the start of the spectrumList
curLine = fgetl(fidRead);

if strcmp(datasetType, 'profile spectrum')
    axeMZ = [];
    count = 1;
    while ~strcmp(LDR.label, '/run') && count <= scanCount
        if strcmp(LDR.label, 'spectrum'), boolArray = 1; end
        if strcmp(LDR.label, 'cvParam') && ~isempty(LDR.attributes)
            [~, ~, field] = fieldfind( LDR.attributes, 'accession');
            
        end
        if strcmp(LDR.label, 'binary')
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
                if options.display
                    disp(['Checking the m/z axe ', num2str(count), ...
                        ' out of ', num2str(scanCount)])
                end
                axeMZ = unique([axeMZ, round(mzValue, 4)]);
            else
                boolArray = 1;
                intValue = output;
                count = count + 1;
            end
        else
        end
        
        
        curLine = fgetl(fidRead);
        if ~ischar(curLine)
            error('myApp:endOfFile', ...
                'The mzML file is not complete')
        end
        [LDR, curLine] = getMZMLCamp(curLine, fidRead);
    end
end



fseek(fidRead, keptPl, 'bof') % Go back at the start of the spectrumList
curLine = fgetl(fidRead);

% 4. EXTRACT DATA FOR EACH SCAN
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
        if length(TICP) == scanCount, break; end
        if axeX(end) >= finneeStc.info.parameters.xMin && ...
                axeX(end) <= finneeStc.info.parameters.xMax
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
                
                MS =  [mzValue' intValue'];
                if strcmp(datasetType, 'profile spectrum')
                    MS_new = axeMZ';
                    MS_new(:,2) = 0;
                    Lia = ismember(axeMZ, round(MS(:,1), 4));
                    assignin('base', 'axeMZ', axeMZ)
                    assignin('base', 'MS', MS)
                    MS_new(Lia, 2) = MS(:,2);
                    MS = MS_new;
                end
                
                % calculates profiles and limits
                if mzMin > min(MS(:,1)), mzMin = min(MS(:,1)); end
                if mzMax < max(MS(:,1)), mzMax = max(MS(:,1)); end
                if intMin > min(MS(:,2)), intMin = min(MS(:,2)); end
                if intMax < max(MS(:,2)), intMax = max(MS(:,2)); end
                TICP(end+1) = sum(MS(:,2));
                [BPP(end+1), indMax] = max(MS(:,2));
                mzBPP(end+1) = MS(indMax, 2);
                % write
                posIni =  ftell(fidWriteDat);
                fwrite(fidWriteDat, ...
                    MS, 'double');
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
warning('on','all')

% 5. CONCLUSIONS AND DATA RECORDING
save2struc()
finneeStc.dataset{1}.info.errors = {};
save(fullfile(finneeStc.info.parameters.folderOut, ...
    [finneeStc.info.parameters.fileID '.mat']), 'finneeStc')
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
                curLine = fgetl(fidRead)
                [LDR, ~] = getMZMLCamp(curLine, fidRead);
                if strcmp(endLabel, LDR.label)
                    break
                end
                if ~isfield(outStrc.(parentLabel), LDR.label)
                    outStrc.(parentLabel).(LDR.label) = {};
                end
                s = addChildNode(LDR);
                
                assignin('base', 's', s)
                outStrc.(parentLabel).(LDR.label){end+1} = s.(LDR.label);
            end
        end
        
    end

% 2. NESTED-FUNCTION:  save2struc()
    function save2struc()
    % Save description of dataset
    finneeStc.dataset{1}.name = 'Original dataset';
    finneeStc.dataset{1}.dateOfCreation = datetime;
    finneeStc.dataset{1}.info = info;
    finneeStc.dataset{1}.info.parameters = ...
        finneeStc.info.parameters;
    finneeStc.dataset{1}.info.errors = {};
    finneeStc.dataset{1}.description.mzStart = mzMin;
    finneeStc.dataset{1}.description.mzEnd = mzMax;
    finneeStc.dataset{1}.description.intMin = intMin;
    finneeStc.dataset{1}.description.intMax = intMax;
    finneeStc.dataset{1}.description.timeStart = axeX(1);
    finneeStc.dataset{1}.description.timeEnd = axeX(length(axeX));
    finneeStc.dataset{1}.description.dataFormat = datasetType;
    finneeStc.dataset{1}.indexInDat = MSIndex;

    finneeStc.dataset{1}.axes.time.values = axeX';
    finneeStc.dataset{1}.axes.time.label = timeLabel;
    finneeStc.dataset{1}.axes.time.unit = timeUnit;
    finneeStc.dataset{1}.axes.mz.values = [];
    finneeStc.dataset{1}.axes.mz.label = mzLabel;
    finneeStc.dataset{1}.axes.mz.unit = mzUnit;
    finneeStc.dataset{1}.axes.intensity.values = [];
    finneeStc.dataset{1}.axes.intensity.label = intLabel;
    finneeStc.dataset{1}.axes.intensity.unit = intUnit;
   
    
    % Record profiles
    % ** TICP
    finneeStc.dataset{1}.trace{1}.name = 'Total Ion Current Profile (dataset 1)';
    finneeStc.dataset{1}.trace{1}.dateOfCreation = datetime;
    finneeStc.dataset{1}.trace{1}.code = 'TIP';
    finneeStc.dataset{1}.trace{1}.plotType = 'profile';
    finneeStc.dataset{1}.trace{1}.axeX.label = timeLabel;
    finneeStc.dataset{1}.trace{1}.axeX.unit = timeUnit;
    finneeStc.dataset{1}.trace{1}.axeY.label = intLabel;
    finneeStc.dataset{1}.trace{1}.axeY.unit = intUnit;
    finneeStc.dataset{1}.trace{1}.indexInDat  = [ftell(fidWriteDat), 0, 2];
    assignin('base', 'axeX', axeX)
    assignin('base', 'TICP', TICP)
    % WITH msconvert and agilent TICP is longer than axeX in the last
    % scan??? Correction add for this case
    % fwrite(fidWriteDat, [axeX' TICP'], 'double');
    indEnd = min(length(axeX), length(TICP));
    fwrite(fidWriteDat, [axeX(1:indEnd)' TICP(1:indEnd)'], 'double');
    finneeStc.dataset{1}.trace{1}.indexInDat(2) = ftell(fidWriteDat);
    
    % ** BPP
    finneeStc.dataset{1}.trace{2}.name = 'Base Peak Profile (dataset 1)';
    finneeStc.dataset{1}.trace{2}.dateOfCreation = datetime;
    finneeStc.dataset{1}.trace{1}.code = 'BPP';
    finneeStc.dataset{1}.trace{2}.plotType = 'profile';
    finneeStc.dataset{1}.trace{2}.axeX.label = timeLabel;
    finneeStc.dataset{1}.trace{2}.axeX.unit = timeUnit;
    finneeStc.dataset{1}.trace{2}.axeY.label = intLabel;
    finneeStc.dataset{1}.trace{2}.axeY.unit = intUnit;
    finneeStc.dataset{1}.trace{2}.indexInDat  = [ftell(fidWriteDat), 0, 2];
    fwrite(fidWriteDat, [axeX(1:indEnd)' BPP(1:indEnd)'], 'double');
    finneeStc.dataset{1}.trace{2}.indexInDat(2) = ftell(fidWriteDat);
    
    % ** mzBPP
    finneeStc.dataset{1}.trace{3}.name = 'm/z @ Base Peak (dataset 1)';
    finneeStc.dataset{1}.trace{3}.dateOfCreation = datetime;
    finneeStc.dataset{1}.trace{1}.code = 'mzBPP';
    finneeStc.dataset{1}.trace{3}.plotType = 'profile';
    finneeStc.dataset{1}.trace{3}.axeX.label = timeLabel;
    finneeStc.dataset{1}.trace{3}.axeX.unit = timeUnit;
    finneeStc.dataset{1}.trace{3}.axeY.label = mzLabel;
    finneeStc.dataset{1}.trace{3}.axeY.unit = mzUnit;
    finneeStc.dataset{1}.trace{3}.indexInDat  = [ftell(fidWriteDat), 0, 2];
    fwrite(fidWriteDat, [axeX(1:indEnd)' mzBPP(1:indEnd)'], 'double');
    finneeStc.dataset{1}.trace{3}.indexInDat(2) = ftell(fidWriteDat);
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
parameters.prec4mz = '%0.4f';
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
            case 'dispOff'
                options.display = 0;
                SFi = SFi +1;
            case 'overwrite'
                options.ifExist = 'Overwrite';
                SFi = SFi + 1;
            case 'precision'
                parameters.prec4mz = ['%0.', ...
                    num2str(vararginIn{SFi+1}(1), 0), 'f'];
                SFi = SFi + 2;
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

