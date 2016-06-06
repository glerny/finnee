function traceOut = getTrace(finneeStc, address, varargin)
%% DESCRIPTION
% 1. INTRODUCTION
% GETTRACE is used to retrieve a trace (profile, MS spectra, ...) that has
% been save in the  finneeStc. dataOut is the mxn array that contain all
% data related to the specified trace.
%
% 2. INPUT PARAMETERS:
%   . required. GETTRACE requires at least 2 parameters
%       finneeStc
%           is the finnee structure 
%       address
%           address is the location of the index in the structure. the
%           format should be: 'trace@dataset' 
%       
%   .optionals. VARARGIN describes the optional paramters.
%       'noFig'     
%           Will not plot the trace
%
% 3. OUPUT PARAMETER
%   traceOut is a strcuture that will contain the data (x and y) as
%   well as the units, labels and title.
%
% 3. EXAMPLES:
%       dataOut = plotTrace(finneeStc, '3@1')
%
% 4. COPYRIGHT
% Copyright 2014-2015 G. Erny (guillaume@fe.up,pt), FEUP, Porto, Portugal

%% CORE OF THE FUNCTION
% 1. INITIALISATION
info.function.functionName = 'getTrace';
info.function.description{1} = 'get the trace recorded in finneeStc';
info.function.matlabVersion = '8.5.0.197613 (R2015a)';
info.function.version = '14/01/2016';
info.function.ownerContact = 'guillaume@fe.up.pt';

[parameters, options] = initFunction(nargin, finneeStc, address, varargin );
%INITFUNCTION used to verify the entries and load the optional and
% complusory parameters

m = parameters.dataset;
n = parameters.trace;

fidReadTra = fopen(finneeStc.path2dat, 'rb');
index = finneeStc.dataset{m}.trace{n}.indexInDat;
fseek(fidReadTra,  index(1), 'bof');
traceOut.title = finneeStc.dataset{m}.trace{n}.name;
traceOut.data = ...
    fread(fidReadTra, [(index(2)-index(1))/(8*index(3)) index(3)], 'double');
traceOut.plotType = finneeStc.dataset{m}.trace{n}.plotType;
traceOut.axes.axeX = finneeStc.dataset{m}.trace{n}.axeX;
traceOut.axes.axeY = finneeStc.dataset{m}.trace{n}.axeY;
traceOut.info = info;

if options.display
    switch traceOut.plotType
        case 'profile'
            plot(traceOut.data(:,1), traceOut.data(:,2));
            title(traceOut.title);
            xlabel([traceOut.axes.axeX.label, ' / ', traceOut.axes.axeX.unit]);
            ylabel([traceOut.axes.axeY.label, ' / ', traceOut.axes.axeY.unit]);
        case 'stem'
            stem(traceOut.data(:,1), traceOut.data(:,2), 'Marker', 'none');
            title(traceOut.title);
            xlabel([traceOut.axes.axeX.label, ' / ', traceOut.axes.axeX.unit]);
            ylabel([traceOut.axes.axeY.label, ' / ', traceOut.axes.axeY.unit]);
        case 'bar'
            bar(traceOut.data(:,1), traceOut.data(:,2));
            title(traceOut.title);
            xlabel([traceOut.axes.axeX.label, ' / ', traceOut.axes.axeX.unit]);
            ylabel([traceOut.axes.axeY.label, ' / ', traceOut.axes.axeY.unit]);
    end
end

%% NESTED FUNCTIONS
end
%% SUB FUNCTIONS
% 1. INITFUNCTION
% Function that get the input argument and check for errors
function [parameters, options] = ...
    initFunction(narginIn, finneeStc, address, vararginIn )

options.display = 1;
% 1.1. Check for obligatory parameters
if narginIn < 2 % check the number of input parameters
    error('myApp:argChk', ...
        ['Wrong number of input arguments. \n', ...
        'Type help plotTrace for more information']);
elseif ~ischar(address)
    error('myApp:argChk', ...
        ['address shoud be a string. \n', ...
        'Type help MSdata2struct for more information']);
elseif ~isstruct(finneeStc)
    error('myApp:argChk', ...
        ['finneeStc shoud be a structure. \n', ...
        'Type help MSdata2struct for more information']);
end

% 1.2. Check for option
if  narginIn > 2
    SFi = 1;
    length(vararginIn)
    while SFi <= length(vararginIn)
        switch vararginIn{SFi}
            case 'noFig'
                options.display = 0;
                SFi = SFi + 1;
            otherwise
                error('myApp:argChk', ...
                    [vararginIn{SFi} ' is not a recognized PropertyName'])
        end
    end
end

% 1.3. Decifer address and check for errors

list = regexp(address, '@','split');
tgtDataset = str2double(list{2});
tgtTrace =  str2double(list{1});
% check for error
if isempty(tgtDataset)
     error('myApp:argChk', ...
         [address ' is not a recognized address'])
elseif tgtDataset > length(finneeStc.dataset)
    error('DIY')
end

if isempty(tgtTrace)
    error('myApp:argChk', ...
        [address ' is not a recognized address'])
elseif tgtTrace > length(finneeStc.dataset{tgtDataset}.trace)
             error('DIY')
end
parameters.dataset = tgtDataset; parameters.trace = tgtTrace;
        
end

