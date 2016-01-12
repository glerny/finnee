function finneeStc = doHACA(finneeStc, dataset, varargin)
%% DESCRIPTION
% 1. INTRODUCTION
% DOHACA allow to build a hierarchical strcutures by classifing every
% individual PIP contained in a 'ionic profile' dataset via their
% correlation coefficient.
%
% 2. PARAMETERS:
%   .required. DOHACA requires at least 2 parameters
%       finneeStc
%           is the finnee structure that contain information about the run
%           and link and indexation of the associated dat file. The
%           strcuture should have been create by function such as 
%           MZML2STRUCT
%       dataset
%           datasetis the indice to the targeted dataset (i.e. in
%           finneeStc.dataset{m}, where m is the target dataset). The
%           dataset should be a 'ionic profile' dataset
%
%   .optionals. VARARGIN describes the optional paramters.  
%       'maxProfiles' followed by an integer (default 2000)
%           Allow to limit the number of PIP based on their maximum
%           intensuty. This allow to speed up the compuational time.
%       'minCorCoef' followed by a number (default 0.7)
%           Stop the HACA when the correlation is lower than 0.7, this
%           mainly for size reason
%       'minPIPperCluster' follwoed by an integer (default 2)
%           Only record clusters with 'minPIPperCluster' PIP
%
% 3. EXAMPLES:
%	finneeStc = doHACA(finneeStc, 2)
%
% 4. COPYRIGHT
% Copyright 2014-2015 G. Erny (guillaume@fe.up.pt), FEUP, Porto, Portugal
%

%% CORE OF THE FUNCTION
% 1. INITIALISATION
info.functionName = 'doHACA';
info.description{1} = 'Calculate the herarchical structure of '' dataset';
info.matlabVersion = '8.5.0.197613 (R2015a)';
info.version = '09/07/2015_gle01';
info.ownerContact = 'guillaume@fe.up.pt';
[parameters, options] = ...
    initFunction(nargin,  finneeStc, dataset,  varargin );
%INITFUNCTION - sub function used to verify the entries and load the optional and
% compulsory and optional parameters

m = parameters.dataset;

% 2. CHECKING THE DATA TYPE
switch finneeStc.dataset{m}.description.dataFormat
    case 'ionic profile'
        fidReadDat = fopen(finneeStc.dataset{m}.description.path2DatFile, 'rb');
        fom.label = finneeStc.dataset{m}.description.fom.label;
        index = finneeStc.dataset{m}.description.fom.data;
        fseek(fidReadDat, index(1), 'bof');
        fom.data = sortrows(fread(fidReadDat, ...
            [(index(2)-index(1))/(index(3)*8), index(3)], 'double'), -5);
        index = finneeStc.dataset{m}.description.axe;
        fseek(fidReadDat, index(1), 'bof');
        axeX = fread(fidReadDat, [(index(2)-index(1))/(index(3)*8), ...
            index(3)], 'double');
        matrixOfPIP = zeros(length(axeX), ...
            min(length(fom.data(:,1)), parameters.maxPIP));
        
        
        % 2. INITIALIZING THE HIERACHICAL STRUCTURE, CALCULATING THE MATRIX
        % OF CORRELATION
        if options.display, disp('Calculating matrix of correlations'); end
        for ii = 1:length(matrixOfPIP(1,:))
            index = finneeStc.dataset{m}.description.index2DotDat(fom.data(ii,1), :);
            fseek(fidReadDat, index(1), 'bof');
            curPIP = fread(fidReadDat, [(index(2)-index(1))/(index(3)*8), ...
                index(3)], 'double');
            % doFilter should not be used at the moment will be tested
            % latter
            if parameters.doFilter
                curPIP(:,2) = smooth(curPIP(:,2),'sgolay',3);
            end
            matrixOfPIP(curPIP(:,1), ii) = curPIP(:,2);
            indOfPIP(ii) = fom.data(ii,1);
        end
        fclose(fidReadDat);
        
        clustPlot = [];
        assignin('base', 'matrixOfPIP', matrixOfPIP)
        assignin('base', 'indOfPIP', indOfPIP)
        data4cluster = {};
        STOPMEHERE
        
        % 3. FILLING THE HIERARCHICAL STRUCTURE
        while 1
            %             if options.display
            
            %             end
            R = round(corrcoef(matrixOfPIP), 2);
            max_CC = max(max(triu(R,1)));
            [linMin, colMin] = find(triu(R,1) == max_CC);
            [~, indL] = min(linMin);
            prof1 = matrixOfPIP(:,linMin(indL));
            if ~isempty(clustPlot)
                R2 = round(corrcoef([prof1, clustPlot]), 2);
                if max(R2(1,2:end))+ 0.01 >= max_CC
                    [~, indconcat] = max(R2(1,2:end));
                    if length(indconcat) > 1
                        NOPODESER
                    end
                    data4cluster{indconcat}.listOfPIP(end+1) = ...
                        indOfPIP(linMin(indL));
                    matrixOfPIP(:,linMin(indL)) = [];
                    indOfPIP(linMin(indL)) = [];
                else
                    clustPlot(:, end+1) = (matrixOfPIP(:,linMin(indL)) + ...
                        matrixOfPIP(:,colMin(indL)));
                    data4cluster{end+1}.listOfPIP = [indOfPIP(linMin(indL)), ...
                        indOfPIP(colMin(indL))];
                    data4cluster{end}.mergCC = max_CC;
                    matrixOfPIP(:,[linMin(indL),colMin(indL)]) = [];
                    indOfPIP([linMin(indL),colMin(indL)]) = [];
                end
            else
                clustPlot(:, end+1) = (matrixOfPIP(:,linMin(indL)) + ...
                    matrixOfPIP(:,colMin(indL)));
                data4cluster{end+1}.listOfPIP = [indOfPIP(linMin(indL)), ...
                    indOfPIP(colMin(indL))];
                data4cluster{end}.mergCC = max_CC;
                matrixOfPIP(:,[linMin(indL),colMin(indL)]) = [];
                indOfPIP([linMin(indL),colMin(indL)]) = [];
            end
            fprintf('Clustering step : free profile:  %d \n\t total clusters: %d \n',...
                length(matrixOfPIP(1,:)), length(clustPlot(1,:)))
        end
end
            
%             % 4.  SAVING THE HIERARCHICAL STRUCTURE
%             if ~isfield(finneeStc.dataset{m}, 'clusterAnalysis')
%                 finneeStc.dataset{m}.clusterAnalysis = {};
%             end
%             finneeStc.dataset{m}.clusterAnalysis{end+1}.infoFunctionUsed...
%                 .info = info;
%             finneeStc.dataset{m}.clusterAnalysis{end}.infoFunctionUsed...
%                 .parameters = parameters ;
%             finneeStc.dataset{m}.clusterAnalysis{end}.infoFunctionUsed...
%                 .error = {};
%             finneeStc.dataset{m}.clusterAnalysis{end}.description...
%                 .dateOfCreation = clock;
%             
%             fidWriteDat = ...
%                 fopen(finneeStc.dataset{m}.description.path2DatFile, 'ab');
%             cl2add.step = {};
%             for ii = 1:length(clusters.step)
%                 if options.display,
%                     fprintf('Saving cluster : %d / %d \n', ii, ...
%                         length(clusters.step))
%                 end
%                 cl2add.step{end+1}.HL = clusters.step{ii}.HL;
%                 cl2add.step{end}.index2DotDat = [ftell(fidWriteDat), 0, 1];
%                 data2write = zeros(1, length(clusters.step));
%                 for jj = 1:length(clusters.step{ii}.cluster)
%                     if length(clusters.step{ii}.cluster{jj}) >= parameters.PIPperCl
%                         data2write(end+1, 1:length(clusters.step{ii}.cluster{jj})) = ...
%                             fom.data(clusters.step{ii}.cluster{jj}, 1);
%                     end
%                 end
%                 data2write(1,:) = [];
%                 ind2rem = sum(data2write) == 0;
%                 data2write(:,ind2rem) = [];
%                 if isempty(data2write)
%                     cl2add.step{end}.index2DotDat(2) = ftell(fidWriteDat);
%                     cl2add.step{end}.index2DotDat(3) = 0;
%                 else
%                     fwrite(fidWriteDat, data2write, 'double');
%                     cl2add.step{end}.index2DotDat(2) = ftell(fidWriteDat);
%                     cl2add.step{end}.index2DotDat(3) = length(data2write(1,:));
%                 end
%             end
%             
%             finneeStc.dataset{m}.clusterAnalysis{end}.HStructure =...
%                 cl2add;
%             fclose(fidWriteDat);
%             save(fullfile(finneeStc.infoFunctionUsed.parameters.folderOut, ...
%                 [finneeStc.infoFunctionUsed.parameters.fileID '.fin']), 'finneeStc', '-mat')
%         end
%     otherwise
%         error('the dataset should be ionic profile')
% end

        %% NESTED FUNCTIONS
end

%% SUB FUNCTIONS
% 1. INITFUNCTION
% Function that get the input arguments and check for errors
    function [parameters, options] = ...
            initFunction(narginIn, finneeStc, dataset, vararginIn )
        
        % defaults parameters
        parameters.minCC = 0.7;
        parameters.maxPIP = 2000;
        parameters.PIPperCl = 2;
        options.display = 1;
        parameters.dataset = dataset;
        parameters.doFilter = 0;
        
        % 1.1. Check for obligatory parameters
        if narginIn < 2 % check the number of input parameters
            error('myApp:argChk', ...
                ['Wrong number of input arguments. \n', ...
                'Type help plotTrace for more information']);
        elseif  ~isnumeric(dataset)
            error('myApp:argChk', ...
        ['DATASET shoud be a string. \n', ...
        'Type help MSdata2struct for more information']);
elseif ~isstruct(finneeStc)
    error('myApp:argChk', ...
        ['finneeStc shoud be a structure. \n', ...
        'Type help MSdata2struct for more information'])
end
% 1.2. Check for optional parameters
if  narginIn > 2
    SFi = 1;
    length(vararginIn)
    while SFi <= length(vararginIn)
        switch vararginIn{SFi}
            case 'maxProfiles'
                parameters.maxPIP = vararginIn{SFi+1};
                SFi = SFi + 2;
            case 'minCorCoef'
                parameters.minCC = vararginIn{SFi+1};
                SFi = SFi + 2;
            case 'minPIPperCluster'
                parameters.minCC = vararginIn{SFi+1};
                parameters.PIPperCl = SFi + 2;
            otherwise
                error('myApp:argChk', ...
                    [vararginIn{SFi} ' is not a recognized PropertyName'])
        end
    end
end
end


