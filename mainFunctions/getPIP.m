function [finneeStc, PIPOut] = getPIP( finneeStc, dataset, intThresh, diffMass, ptsPerPeak, varargin)
%% DESCRIPTION
% 1. INTRODUCTION
%   GETPIP allows to extract from a 'centroid spectrum' dataset pure ion profiles
%   (PIP) using, as input parameters, an intensity threshold (intThresh) allowing to 
%   separate noise from relevant information and the maximum mass difference between
%   datapoints in successive scans to potentially belong to the same ionic 
%   profile (diffMass). 
% 2. INPUT PARAMETERS
%   .required. GETPIP requires at least 5 parameters
%       finneeStc
%           is the structure that contains information about the run
%           as well as link and indexation to the associated dat file. 
%       dataset
%           dataset is the indice to the targeted dataset (i.e.
%           finneeStc.dataset{m}, where m is the target dataset). The
%           dataset should be a 'centroid spectrum' dataset.
%       intThresh
%           Only datapoints of intensity higher than intThresh will be
%           consider in a first instance to find the mass range of interest 
%           (MRI). MRI are determined by projecting all datapoints whose intensity
%           are higher than intThresh in a m/z axe. MRI are defined as any succession
%           of points whose m/z values does not differ by more that a set value (diffMass) 
%           to its closest neighbourg. 
%       diffMass
%           Maximum difference between two succesive datapoints in the m/z projection
%           to be potentially part of the same ionic profile.
%       ptsPerPeak
%           Each MRI are then treated to regroup any points whose scan number are
%           concomitant. Sequences of more than ptsPerPeak are assumed to be part of a 
%           pure ion profile (PIP).
%
%   .optionals. VARARGIN describes the optional parameters.  
%       'mzInt' followed by a 2x1 array
%           Allow to only consider the datapoints with m/z are within the
%           defined interval.
%       'timeInt' followed by a 2x1 array
%           Allow to only consider the datapoints with time (scan number)
%           are within the defined interval.
%
% 3. OUTPUT PARAMETERS:
%
% 4. EXAMPLES:
%   [finneeStc, PIPOut] = getPIP( finneeStc, 2, 25, 0.005, 4)
%
% 5. COPYRIGHT
%   Copyright 2015-2016 G. Erny (guillaume@fe.up.pt), FEUP, Porto, Portugal
%

%% CORE OF THE FUNCTION
% 1. INITIALISATION
info.function.functionName = 'getPIP';
info.function.description{1} = 'get the PIP from the original ''centroid spectrum'' dataset';
info.function.matlabVersion = '8.5.0.197613 (R2015a)';
info.function.version = '18/01/2016';
info.function.ownerContact = 'guillaume@fe.up.pt';

[parameters, options] = ...
    initFunction(nargin,  finneeStc, dataset, intThresh, diffMass, ...
    ptsPerPeak,  varargin );
%INITFUNCTION used to verify the entries and load the optional and
% compulsory parameters

m = parameters.dataset;
iT = parameters.intThresh;
d = parameters.diffMass;
p = parameters.ptsPerPeak;


fidReadDat = fopen(finneeStc.path2dat, 'rb');

% 2. FINDING SEQUENCES OF POINT THAT MAY BELONG TO A PIP
% 2.1. Filtering out dara whith I < intThres.
%             sumMS is a 3xn array that will contain
%             all remaining n points with in first column 
%             the scan number, the second the intensity 
%             and the third the m/z values
switch finneeStc.dataset{m}.description.dataFormat
    case 'centroid spectrum'    % !THIS FUNCTION ONLY WORK WITH CENTROID 
                                % SPECTRUM DATASET!
                                
        axeX = finneeStc.dataset{m}.axes.time.values;
        indTimeStt = findCloser(parameters.xMin, axeX);
        indTimeEnd = findCloser(parameters.xMax, axeX);
        sumMS = []; 
                    
        for ii = indTimeStt:indTimeEnd % for each MS scan do...
            if options.display 
                fprintf('processing scan %d out of %d scans \n', ii, length(axeX))
            end
            index = finneeStc.dataset{m}.indexInDat(ii, :);
            fseek(fidReadDat, index(1), 'bof');
            MS = ...
                fread(fidReadDat, [(index(2)-index(1))/(index(3)*8), index(3)], 'double');
				% MS is ths centroid MS scan ii with in first column the m/z values and in
				% second columns the corresponding intensity counts
				
            ind2rem = MS(:,1) < parameters.mzMin |...
                MS(:,1) > parameters.mzMax;
            MS(ind2rem, :) = []; 
			% Step needed to remove the pts ourside the m/z zone of interest
			
            ind2rem = MS(:,2) < iT;
            MS(ind2rem, :) = [];
			% filtering to remove data below intThresh
			
	        MS(:,3) = ii;
			% adding scan number
			
            if ~isempty(MS)
                if isempty(sumMS)
                    sumMS = [MS(:,3) MS(:,2) MS(:,1)];
                else
                    sumMS = [sumMS; [MS(:,3) MS(:,2) MS(:,1)]];
                end
            end
			% concatening new MS data to sumMS
			
        end
    otherwise
        fclose(fidReadDat);
        error('dataset should be ''centroid spectrum'' dataset')
end
fclose(fidReadDat);

% 2.2. Finding sequences of pts that may belong to PIP
sumMS = sortrows(sumMS, 3);
% sorting sumMS by m/z values.

diffMZ = diff(sumMS(:,3));
ind2cutMZ = [0; find(diffMZ > d); length(sumMS(:,1))];
% ind2cutMZ contains all indice ind such as  
% sumMS(ind+1, 3) - sumMS(ind, 3) > diffMass
 
ind2cut2MZ = find(diff(ind2cutMZ) >= p);
% ind2cut2MZ verify that the previous intervals
% contain more than ptsPerPeak

MII = {};

for ii = 1:length(ind2cut2MZ)
   provData = sortrows...
       (sumMS(ind2cutMZ(ind2cut2MZ(ii))+1:ind2cutMZ(ind2cut2MZ(ii)+1),:), 1);
   % provData contains all data points related to the mass range of interest (MRI) ii. 
   % Those are series of at least ptsPerPeak datapoints, each of them has its m/z that 
   % does not differ by more than diffMass from its closest neighbourgh.
   
   diffTime =  diff(provData(:,1));
   ind2cutT = [0; find(diffTime > 1); length(provData(:,1))];
   ind2cut2T = find(diff(ind2cutT) >= p);
   % ind2cut2T contains the index of series of pts in provData that follow eah other in
   % their scan numbers
   
   for jj = 1:length(ind2cut2T)
       MII{end+1} = ...
           provData(ind2cutT(ind2cut2T(jj))+1:ind2cutT(ind2cut2T(jj)+1),:);
   end
end


% 3.3. Each partial PIP are verify for time doublons
%      (i.e. two points in the MII with the same scan number)
ii = 1;
while 1
    if options.display
        fprintf('Checking for doublons %d out of %d \n', ii, length(MII))
    end
    curData = MII{ii};
    diffTime =  diff(curData(:,1));
    ind = find(diffTime == 0);
    if ~isempty(ind)
        [~, indsq] = max(abs(curData(ind, 3) - curData(ind+1, 3)));
        % If more than one doublon find the doublon with the highest m/z
        % difference
        MII(ii) = [];
        cut1 =  curData(ind(indsq), :); % start two lines
        cut2 =  curData(ind(indsq)+1, :);
        scan2cut = curData(ind(indsq)+1, 1);
        curData(ind(indsq)+1, :) = [];
        curData(ind(indsq) , :) = [];
        
        % Order curData by distance to ind(indsq)
        dist2cut = abs(curData(:,1) - scan2cut);
        [~,I] = sort(dist2cut);
        curData = curData(I,:);
        for jj = 1:length(curData(:,1)) % for each pts put in the lines
            % with the closest m/z
            if abs(curData(jj,3) - mean(cut1(:,3))) > ...
                    abs(curData(jj,3) - mean(cut2(:,3)))
                cut2 = [cut2; curData(jj,:)];
            else
                cut1 = [cut1; curData(jj,:)];
            end
        end
        
        % checking for sequence in cut1 & 2
        cut1 = sortrows(cut1, 1);
        ind2cutT = [0; find( diff(cut1(:,1)) > 1); length( cut1(:,1))];
        ind2cut2T = find(diff(ind2cutT) >= p);
        
        for jj = 1:length(ind2cut2T)
            MII{end+1} = ...
                cut1(ind2cutT(ind2cut2T(jj))+1:ind2cutT(ind2cut2T(jj)+1),:);
        end
        
        cut2 = sortrows(cut2, 1);
        ind2cutT = [0; find( diff(cut2(:,1)) > 1); length( cut2(:,1))];
        ind2cut2T = find(diff(ind2cutT) >= p);
        for jj = 1:length(ind2cut2T)
            MII{end+1} = ...
                cut2(ind2cutT(ind2cut2T(jj))+1:ind2cutT(ind2cut2T(jj)+1),:);
        end
        
    else
        ii = ii + 1;
    end
    if ii > length(MII), break; end
end

% %%%%%%%%%%%%%% ADDING THE 08-10-2015 FIND ionic distribution that belong
% %%%%%%%%%%%%%% to the same ion (t-test)

if parameters.concat.do
    fom.data = zeros(length(MII), 4);
    for ii = 1:length(MII)
        curData = MII{ii};
        % 	% 4.2. Calculating figures of merits
        fom.data(ii, 1) = ii;
        fom.data(ii, 2) = mean(curData(:,3));
        fom.data(ii, 3) = std(curData(:,3));
        fom.data(ii, 4) = length(curData(:,1));
    end
    
    ll = length(fom.data(:,1));
    ResCorr = zeros(ll,1);
    for ii = 1:ll
        if options.display
            fprintf('t-test and concatenated: %d out of %d \n', ii, ll)
        end
        s = ((fom.data(ii,4)-1)*fom.data(ii,3)^2+ (fom.data(:,4)-1).*...
            fom.data(:,3).^2)./(fom.data(ii,4)+fom.data(:,4)-2);
        t = abs(fom.data(ii,2)-fom.data(:,2))./(sqrt(s).*...
            (sqrt(1/fom.data(ii,4)+1/fom.data(:,4)))');
        degFre = fom.data(ii,4) + fom.data(:,4) - 2;
        tcoeff = zeros(ll,1);
        for jj = length(parameters.concat.tdis):-1:1
            indm = find(degFre <= parameters.concat.tdis(jj, 1));
            tcoeff(indm) = parameters.concat.tdis(jj, 2);
        end
        ind = find(t <= tcoeff);
        
        if length(ind) > 1
            ind2rem = find(ind == ii);
            ind(ind2rem) = [];
            ResCorr(ii,1) = ii;
            for jj = 1:length(ind)
                ResCorr(ii, jj+1) = ind(jj);
            end
        end
    end
    
    indpur =  find(ResCorr(:,1) == 0); % pure trace not to mingle
    ResCorr(indpur,:) = [];
    toDell = [];
    
    list = {};
    while ~isempty(ResCorr)
        list1 = [];
        list2 = nonzeros(ResCorr(end,:));
        while length(list1) ~= length(list2)
            list1 = list2;
            list2 = [];
            for ii = 1:length(list1)
                ind2list = find(ResCorr(:,1) == list1(ii));
                list2 = [list2; nonzeros(ResCorr(ind2list,:))];
            end
            list2 = unique(list2);
        end
        list{end+1} = list2;
        ind2rem = [];
        for ii = 1:length(list2)
            ind2test = find(ResCorr(:,1) == list1(ii));
            if ~isempty(ind2test)
                ind2rem(end+1) = ind2test;
            end
        end
        ResCorr(ind2rem, :) = [];
    end


    for ii = 1:length(list)
        newSerie = [];
        for jj = 1:length(list{ii})
            newSerie = [newSerie; MII{list{ii}(jj)}];
        end
        newSerie = sortrows(newSerie, 1);
        toDell = [toDell; list{ii}];
        ind2cut = ...
            [0; find(diff(newSerie(:,1)) > 1+ parameters.concat.snul); ...
            length(newSerie(:,1))];
        for jj = 1:length(ind2cut)-1
            MII{end+1} = newSerie(ind2cut(jj)+1:ind2cut(jj+1), :);
        end
    end
    MII(toDell) = [];
end


% 4. FINDING EDGES AND COMPLETING PARTIAL PIP TO PIP 
	% 4.1. Completing partial PIP
    fom.label ={'PIP', 'meanMass', 'stdMass', 'nbrPts', 'Imax', ...
		'timeAtImax', 'Area', 'center', 'varInt', 'Noise', 'Signal'};
	fom.data = zeros(length(MII), 11);	
for ii = 1:length(MII)
    if options.display
        fprintf('Completing PIP %d out of %d\n', ii, length(MII))
    end
    curData = MII{ii};
	
    % getting tailing points
    while 1
        if curData(1, 1) - 1 <1, break; end
        MZ1 = min(curData(:, 3)) - d/2;
        MZ2 = max(curData(:, 3)) + d/2;
        spectraOut = getSpectra(finneeStc, m, axeX(curData(1, 1) - 1), ...
            'mzInt', [MZ1 MZ2], 'noFig');
        if isempty(spectraOut.data), break; end
        [~, d2k] = min(abs(spectraOut.data(:,1)-  mean(curData(:, 3))));
        if spectraOut.data(d2k, 2) > curData(1,2), break; end
        curData = [[curData(1, 1) - 1, spectraOut.data(d2k, 2),...
            spectraOut.data(d2k, 1)]; curData];
    end
        
	% Getting fronting points
    while 1
        if curData(end, 1) + 1 > length(axeX), break; end
        MZ1 = min(curData(:, 3)) - d/2;
        MZ2 = max(curData(:, 3)) + d/2;
        spectraOut = getSpectra(finneeStc, m, axeX(curData(end, 1) + 1), ...
            'mzInt', [MZ1 MZ2], 'noFig');
        if isempty(spectraOut.data), break; end
        [~, d2k] = min(abs(spectraOut.data(:,1)-  mean(curData(:, 3))));
        if spectraOut.data(d2k, 2) > curData(end,2), break; end
        curData = [curData; ...
            [curData(end, 1) + 1, spectraOut.data(d2k, 2), ...
            spectraOut.data(d2k, 1)]];
    end
    MII{ii} = curData;
	
	% 4.2. Calculating figures of merits
    fom.data(ii, 1) = ii;
    fom.data(ii, 2) = mean(curData(:,3));
    fom.data(ii, 3) = std(curData(:,3));
    fom.data(ii, 4) = length(curData(:,1));
    [fom.data(ii, 5), indIMax] = max(curData(:,2));
    fom.data(ii, 6) = axeX(curData(indIMax, 1));
    fom.data(ii, 7) =  trapz(axeX(curData(:, 1)), curData(:, 2));
    fom.data(ii, 8) = trapz(axeX(curData(:, 1)), ...
        axeX(curData(:, 1)).*curData(:, 2))/fom.data(ii, 7);
    fom.data(ii, 9) =  var(curData(:,2));
    Ibckg = curData(:,2) <= 3*intThresh;
    fom.data(ii, 10) = std(curData(Ibckg,2));
    fom.data(ii, 11) = fom.data(ii, 5)-mean(curData([1,2, end-1, end],2));    
end

PIPOut.title = ['''ionic profile'' dataset of the ''centroid spectrum'' dataset ', ...
            num2str(m)];
PIPOut.profiles = MII;
PIPOut.FOM.label = fom.label;
PIPOut.FOM.data = fom.data;
PIPOut.axes = finneeStc.dataset{m}.axes;

saveAndPlot()
save(fullfile(finneeStc.info.parameters.folderOut, ...
    [finneeStc.info.parameters.fileID '.mat']), 'finneeStc')

%% NESTED FUNCTIONS
% 1. SAVEANDPLOT
    function saveAndPlot()
        finneeStc.dataset{end+1}.name = PIPOut.title;
        finneeStc.dataset{end}.dateOfCreation = datetime;
        finneeStc.dataset{end}.info = info;
        finneeStc.dataset{end}.info.parameters = parameters;
        finneeStc.dataset{end}.info.errors = {};

        timeLabel = finneeStc.dataset{m}.axes.time.label;
        timeUnit = finneeStc.dataset{m}.axes.time.unit;
        mzLabel = finneeStc.dataset{m}.axes.mz.label;
        mzUnit = finneeStc.dataset{m}.axes.mz.unit;
        intLabel = finneeStc.dataset{m}.axes.intensity.label;
        intUnit = finneeStc.dataset{m}.axes.intensity.unit;
        
        
        % info dataset
        finneeStc.dataset{end}.description.mzStart = 0;
        finneeStc.dataset{end}.description.mzEnd = inf;
        finneeStc.dataset{end}.description.intMin = 0;
        finneeStc.dataset{end}.description.intMax = inf;
        finneeStc.dataset{end}.description.timeStart = 0;
        finneeStc.dataset{end}.description.timeEnd = inf;
        finneeStc.dataset{end}.description.dataFormat = 'ionic profile';
        finneeStc.dataset{end}.indexInDat = [];
        finneeStc.dataset{end}.axes = finneeStc.dataset{m}.axes;
        finneeStc.dataset{end}.trace = {};
        
        % saving new axe, fom, PIP and calculating profiles
        fidWriteDat = fopen(finneeStc.path2dat, 'ab');
        mzStart = inf; mzEnd = 0; intMin = inf; intMax = 0;
        timeStart = inf; timeEnd = 0;
        
       
        [TICP, BPP, mzBPP] = deal(axeX);
        TICP(:,2) = 0; BPP(:,2) = 0; mzBPP(:,2) = 0;
        finneeStc.dataset{end}.FOM.label = fom.label;
        finneeStc.dataset{end}.FOM.data = fom.data;
        
        for NF2 = 1:length(MII)
            curData = MII{NF2};
            TICP(curData(:,1), 2) = TICP(curData(:,1), 2) + curData(:,2);
            ind2add = BPP(curData(:,1), 2) < curData(:,2);
            BPP(curData(ind2add,1), 2) = curData(ind2add,2);
            mzBPP(curData(ind2add,1), 2) = curData(ind2add,3);
            finneeStc.dataset{end}.indexInDat(NF2,1) = ftell(fidWriteDat);
            fwrite(fidWriteDat, curData, 'double');
            finneeStc.dataset{end}.indexInDat(NF2,2) = ftell(fidWriteDat);
            finneeStc.dataset{end}.indexInDat(NF2,3) = length(curData(1,:));
            if mzStart > min(curData(:,3)), mzStart = min(curData(:,3)); end
            if mzEnd < max(curData(:,3)), mzEnd = max(curData(:,3)); end
            if intMin > min(curData(:,2)), intMin = min(curData(:,2)); end
            if intMax < max(curData(:,2)), intMax = max(curData(:,2)); end
            if timeStart > axeX(min(curData(:,1))), timeStart = axeX(min(curData(:,1))); end
            if timeEnd < axeX(max(curData(:,1))), timeEnd = axeX(max(curData(:,1))); end
        end
        finneeStc.dataset{end}.description.mzStart = mzStart;
        finneeStc.dataset{end}.description.mzEnd = mzEnd;
        finneeStc.dataset{end}.description.intMin = intMin;
        finneeStc.dataset{end}.description.intMax = intMax;
        finneeStc.dataset{end}.description.timeStart = timeStart;
        finneeStc.dataset{end}.description.timeEnd = timeEnd;
        
        % saving profiles
        % ** TICP
        finneeStc.dataset{end}.trace{1}.name = ...
            ['Total Ion Current Profile (dataset ', ...
            num2str(length(finneeStc.dataset)), ')'];
        finneeStc.dataset{end}.trace{1}.dateOfCreation = datetime;
        finneeStc.dataset{end}.trace{1}.code = 'TIP';
        finneeStc.dataset{end}.trace{1}.plotType = 'profile';
        finneeStc.dataset{end}.trace{1}.axeX.label = timeLabel;
        finneeStc.dataset{end}.trace{1}.axeX.unit = timeUnit;
        finneeStc.dataset{end}.trace{1}.axeY.label = intLabel;
        finneeStc.dataset{end}.trace{1}.axeY.unit = intUnit;
        finneeStc.dataset{end}.trace{1}.indexInDat  = [ftell(fidWriteDat), 0, 2];
        fwrite(fidWriteDat, TICP, 'double');
        finneeStc.dataset{end}.trace{1}.indexInDat(2) = ftell(fidWriteDat);
        
        % ** BPP
        finneeStc.dataset{end}.trace{2}.name = ...
            ['Base Peak Profile (dataset ', ...
            num2str(length(finneeStc.dataset)), ')'];
        finneeStc.dataset{end}.trace{2}.dateOfCreation = datetime;
        finneeStc.dataset{end}.trace{1}.code = 'BPP';
        finneeStc.dataset{end}.trace{2}.plotType = 'profile';
        finneeStc.dataset{end}.trace{2}.axeX.label = timeLabel;
        finneeStc.dataset{end}.trace{2}.axeX.unit = timeUnit;
        finneeStc.dataset{end}.trace{2}.axeY.label = intLabel;
        finneeStc.dataset{end}.trace{2}.axeY.unit = intUnit;
        finneeStc.dataset{end}.trace{2}.indexInDat  = [ftell(fidWriteDat), 0, 2];
        fwrite(fidWriteDat, BPP, 'double');
        finneeStc.dataset{end}.trace{2}.indexInDat(2) = ftell(fidWriteDat);

        % ** mzBPP
         finneeStc.dataset{end}.trace{3}.name = ...
            ['m/z @ Base Peak (dataset ', ...
            num2str(length(finneeStc.dataset)), ')'];
        finneeStc.dataset{end}.trace{3}.dateOfCreation = datetime;
        finneeStc.dataset{end}.trace{1}.code = 'mzBPP';
        finneeStc.dataset{end}.trace{3}.plotType = 'profile';
        finneeStc.dataset{end}.trace{3}.axeX.label = timeLabel;
        finneeStc.dataset{end}.trace{3}.axeX.unit = timeUnit;
        finneeStc.dataset{end}.trace{3}.axeY.label = mzLabel;
        finneeStc.dataset{end}.trace{3}.axeY.unit = mzUnit;
        finneeStc.dataset{end}.trace{3}.indexInDat  = [ftell(fidWriteDat), 0, 2];
        fwrite(fidWriteDat, mzBPP, 'double');
        finneeStc.dataset{end}.trace{3}.indexInDat(2) = ftell(fidWriteDat);
        
        fclose(fidWriteDat);
        
        origSize =  finneeStc.dataset{m}.indexInDat(end,2) - ...
            finneeStc.dataset{m}.indexInDat(1,1);
        finSize =  finneeStc.dataset{end}.indexInDat(end,2) - ...
            finneeStc.dataset{end}.indexInDat(1,1);
        
       
        TI_ori = getTrace(finneeStc,  ['1@', num2str(m)], 'noFig');
        BP_ori = getTrace(finneeStc,  ['2@', num2str(m)], 'noFig');
        
        figure
        subplot(2,2,1)
        plot(TI_ori.data(:,1), TI_ori.data(:,2));
        hold on
        plot(TICP(:,1), TICP(:,2))
        hold off
        title('Comparaison total ion current profiles')
        xlabel([TI_ori.axes.axeX.label, ' / ', TI_ori.axes.axeX.unit]);
        ylabel([TI_ori.axes.axeY.label, ' / ', TI_ori.axes.axeY.unit]);
        
        subplot(2,2,3)
        plot(TI_ori.data(:,1), TI_ori.data(:,2)-TICP(:,2));
        title('Residual')
        xlabel([TI_ori.axes.axeX.label, ' / ', TI_ori.axes.axeX.unit]);
        ylabel([TI_ori.axes.axeY.label, ' / ', TI_ori.axes.axeY.unit]);
        
        subplot(2,2,2)
        plot(BP_ori.data(:,1), BP_ori.data(:,2));
        hold on
        plot(BPP(:,1), BPP(:,2))
        hold off
        title('Comparaison base peak profiles')
        xlabel([BP_ori.axes.axeX.label, ' / ', BP_ori.axes.axeX.unit]);
        ylabel([BP_ori.axes.axeY.label, ' / ', BP_ori.axes.axeY.unit]);
        
        subplot(2,2,4)
        plot(BP_ori.data(:,1), BP_ori.data(:,2)-BPP(:,2));
        title('Residual')
        xlabel([BP_ori.axes.axeX.label, ' / ', BP_ori.axes.axeX.unit]);
        ylabel([BP_ori.axes.axeY.label, ' / ', BP_ori.axes.axeY.unit]);

        
        fprintf('Summary of the conditions used: \n\n')
        fprintf('\t Intensity threshold: \t \t \t \t \t \t %d \n', intThresh)
        fprintf('\t Maximum m/z difference with the same PIP: \t %d \n', diffMass)
        fprintf('\t Minimum successive points for a PIP: \t \t %d \n\n', ptsPerPeak)
      
        fprintf('Number of PIP extracted from dataset %d:\t %d \n', m, NF2)
        fprintf('Size of the original dataset: \t \t \t %d \n', origSize)
        fprintf('Size of the ''ionic profile'' dataset:\t %d\n', finSize)
        fprintf('Data reduced by: \t \t \t \t \t \t %.1f \n', (origSize-finSize)/origSize*100)
    end
end

%% SUB FUNCTIONS
% 1. INITFUNCTION
% Function that get the input argument and check for errors
function [parameters, options] = ...
    initFunction(narginIn, finneeStc, dataset, intThresh, diffMass, ptsPerPeak, vararginIn )

% 1.1. Check for obligatory parameters
if narginIn < 5 % check the number of input parameters
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

options.display = 1 ;
parameters.intThresh = intThresh;
parameters.diffMass = diffMass;
parameters.ptsPerPeak = ptsPerPeak;
parameters.dataset = dataset;
parameters.mzMin = 0;
parameters.mzMax = inf;
parameters.xMin = 0;
parameters.xMax = inf;
parameters.test = 0;
options.text = 1;
options.save = 1;
parameters.concat.do = 0; %SET TO ZERO THIS PART IS BUGED DOUBLED SOME PROFILES
parameters.concat.tdis = [1, 12.71; 2, 4.30; 3, 3.18; 4, 2.78; 5, 2.57; ...
    6, 2.45; 7, 2.36; 8, 2.31; 9, 2.26; 10, 2.23; 12, 2.18; 14, 2.14; ...
    16, 2.12; 18, 2.10; 20, 2.09; 30, 2.04; 50, 2.01];
parameters.concat.snul = 2;


% 1.2. Check for optional parameters
if  narginIn > 5
    SFi = 1;
    length(vararginIn)
    while SFi <= length(vararginIn)
        switch vararginIn{SFi}
            case 'mzInt'
                if length(vararginIn{SFi+1}) == 2
                    mzInt = sort(vararginIn{SFi+1});
                    parameters.mzMin = mzInt(1);
                    parameters.mzMax = mzInt(2);
                end
                SFi = SFi + 2;
            case 'timeInt'
                if length(vararginIn{SFi+1}) == 2
                    timeInt = sort(vararginIn{SFi+1});
                    parameters.xMin = timeInt(1);
                    parameters.xMax = timeInt(2);
                end
                SFi = SFi + 2;
            case 'noFig'
                options.display  = 0;
                SFi = SFi + 1;
            otherwise
                error('myApp:argChk', ...
                    [vararginIn{SFi} ' is not a recognized PropertyName'])
        end
    end
end
        
end
