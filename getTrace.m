function dataOut = getTrace(finneeStc, adress, varargin)
%% DESCRIPTION
% 1. INTRODUCTION
% GETTRACE is used to retrieve a trace (profile, MS spectra, ...) that has
% been save in the  finneeStc. dataOut is the mxn array that contain all
% data related to the specified trace.
%
% 2. PARAMETERS:
%   . required. GETTRACE requires at least 2 parameters
%       finneeStc
%           is the finnee structure 
%       adress
%           adress is the location of the index in the structure. the
%           format should be: 'trace@dataset' 
%       
%   .optionals. VARARGIN describes the optional paramters.
%       'newFig'    
%           default parameter. Plot the trace in a new figure
%       'inFig' followed by an handle, 
%           will add the figure to an existing figure with the 
%           corresponding handle.
%       'noFig'     
%           Will not plot the trace
%
% 3. EXAMPLES:
%       dataOut = plotTrace(finneeStc, '3@1')
%
% 4. COPYRIGHT
% Copyright 2014-2015 G. Erny (guillaume@fe.up,pt), FEUP, Porto, Portugal

%% CORE OF THE FUNCTION
% 1. INITIALISATION
info.functionName = 'getTrace';
info.description{1} = 'get the trace recorded in finneeStc';
info.matlabVersion = '8.5.0.197613 (R2015a)';
info.version = '03/07/2015_gle01';
info.ownerContact = 'guillaume@fe.up,pt';

[parameters, options] = initFunction(nargin, finneeStc, adress, varargin );
%INITFUNCTION used to verify the entries and load the optional and
% complusory parameters

m = parameters.dataset;
n = parameters.trace;

fidReadTra = fopen(finneeStc.dataset{m}.description.path2DatFile, 'rb');
index = finneeStc.dataset{m}.trace{n}.index2DotDat;
fseek(fidReadTra,  index(1), 'bof');
dataOut = fread(fidReadTra, [(index(2)-index(1))/(8*index(3)) index(3)], 'double');

switch options.display.in
    case 'noFig'
        return
    case 'newFig'
    case 'inFig'
        h = options.display.handle;
end
switch finneeStc.dataset{m}.trace{n}.description.plotType
    case 'profile'
        plot(dataOut(:,1), dataOut(:,2));
        title(finneeStc.dataset{m}.trace{n}.description.name);
        xlabel([finneeStc.dataset{m}.trace{n}.description.axeX.label,...
            ' / ',finneeStc.dataset{m}.trace{n}.description.axeX.unit]);
        ylabel([finneeStc.dataset{m}.trace{n}.description.axeY.label,...
            ' / ',finneeStc.dataset{m}.trace{n}.description.axeY.unit]);
    case 'stem'
        stem(dataOut(:,1), dataOut(:,2), 'Marker', 'none');
        title(finneeStc.dataset{m}.trace{n}.description.name);
        xlabel([finneeStc.dataset{m}.trace{n}.description.axeX.label,...
            ' / ',finneeStc.dataset{m}.trace{n}.description.axeX.unit]);
        ylabel([finneeStc.dataset{m}.trace{n}.description.axeY.label,...
            ' / ',finneeStc.dataset{m}.trace{n}.description.axeY.unit]);
    case 'bar'
        bar(dataOut(:,1), dataOut(:,2));
        title(finneeStc.dataset{m}.trace{n}.description.name);
        xlabel([finneeStc.dataset{m}.trace{n}.description.axeX.label,...
            ' / ',finneeStc.dataset{m}.trace{n}.description.axeX.unit]);
        ylabel([finneeStc.dataset{m}.trace{n}.description.axeY.label,...
            ' / ',finneeStc.dataset{m}.trace{n}.description.axeY.unit]);
end

%% NESTED FUNCTIONS
end
%% SUB FUNCTIONS
% 1. INITFUNCTION
% Function that get the input argument and check for errors
function [parameters, options] = ...
    initFunction(narginIn, finneeStc, adress, vararginIn )

options.display.in = 'newFig';
% 1.1. Check for obligatory parameters
if narginIn < 2 % check the number of input parameters
    error('myApp:argChk', ...
        ['Wrong number of input arguments. \n', ...
        'Type help plotTrace for more information']);
elseif ~ischar(adress)
    error('myApp:argChk', ...
        ['ADRESS shoud be a string. \n', ...
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
            case 'newFig'
                options.display.in = 'newFig';
                SFi = SFi + 1;
            case 'inFig'
                options.display.in = 'inFig';
                options.display.handle = vararginIn{SFi+1};
                SFi = SFi +2;
            case 'noFig'
                options.display.in = 'noFig';
                SFi = SFi + 1;
            otherwise
                error('myApp:argChk', ...
                    [vararginIn{SFi} ' is not a recognized PropertyName'])
        end
    end
end

% 1.3. Decifer adress and check for errors
list = strsplit(adress, '@');
tgtDataset = str2double(list{2});
tgtTrace =  str2double(list{1});
% check for error
if isempty(tgtDataset)
     error('myApp:argChk', ...
         [adress ' is not a recognized adress'])
elseif tgtDataset > length(finneeStc.dataset)
    error('DIY')
end

if isempty(tgtTrace)
    error('myApp:argChk', ...
        [adress ' is not a recognized adress'])
elseif tgtTrace > length(finneeStc.dataset{tgtDataset}.trace)
             error('DIY')
end
parameters.dataset = tgtDataset; parameters.trace = tgtTrace;
        
end

