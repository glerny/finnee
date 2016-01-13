function tellFinnee(finneeStc)
%% DESCRIPTION
% 1. INTRODUCTION
% TELLFINNEE gives informations about the input finneeStc.
%
% 2. PARAMETERS:
%   .required. TELLFINNEE requires at least 1 parameter
%       finneeStc
%           is the finnee structure that contain information about the run
%           and link and indexation of the associated dat file. The
%           strcuture should have been create by function such as 
%           DOMZML2STRUCT
%
% 3. EXAMPLES:
%	tellFinnee(finneeStc)
%
% 4. COPYRIGHT
% Copyright 2015-2016 G. Erny (guillaume@fe.up.pt), FEUP, Porto, Portugal
%

%% CORE OF THE FUNCTION
% 1. General information
%
fprintf('\n\n______________________________________________________________________________\n')
fprintf('\nName: \t\t%s \nFolder: \t%s \n\n', ...
    finneeStc.infoFunctionUsed.parameters.fileID, ...
    finneeStc.infoFunctionUsed.parameters.folderOut);
% 2. File description
fprintf('Original File: %s \n\t(%s  with %s in %s) \n',...
    finneeStc.infoRun.fileDescription.sourceFileList{1}.sourceFile{1}.attributes{2}.field{1},...
    finneeStc.infoRun.fileDescription.sourceFileList{1}.sourceFile{1}.cvParam{2}.attributes{3}.field{1},...
    finneeStc.infoRun.fileDescription.fileContent{1}.cvParam{1}.attributes{3}.field{1},...
    finneeStc.infoRun.fileDescription.fileContent{1}.cvParam{2}.attributes{3}.field{1});
fprintf('Contact Name: \t%s\nContact Organisation: \t%s\n', ...
    finneeStc.infoRun.fileDescription.contact{1}.cvParam{1}.attributes{4}.field{1},...
    finneeStc.infoRun.fileDescription.contact{1}.cvParam{2}.attributes{4}.field{1});
fprintf('Date of recording: \t%s\n\n',...
    finneeStc.infoRun.run.attributes{1, 3}.field{1});
fprintf('Sample Name: \t%s\nSample Description: \n', ...
    finneeStc.infoRun.sampleList.sample{1}.attributes{2}.field{1})
fprintf('\t\t%s\n', ...
    finneeStc.infoRun.sampleList.sample{1}.userParam{1}.attributes{2}.field{:})

% 2. Dataset and traces.
nbrDTS = length(finneeStc.dataset);
fprintf('\nNumber of dataset: \t%d\n', nbrDTS);
for ii = 1:nbrDTS
    fprintf('\t%d.  Name\n', ii)
    fprintf('\t\tData format:\t%s\n\t\tDate of creation:\t%s\n',...
        finneeStc.dataset{ii}.description.dataFormat, ...
        datestr(datetime));
    nbrTRC = length(finneeStc.dataset{ii}.trace);
    fprintf('\t\tNumber of traces: \t%d\n\n', nbrTRC);
    
    for jj = 1:nbrTRC
        fprintf('\t\t%d.%d. %s\n', ii, jj, ...
            finneeStc.dataset{ii}.trace{jj}.description.name);
        fprintf('\t\t\tType:\t%s;\tDate of creation:\t%s\n',...
            finneeStc.dataset{ii}.trace{jj}.description.plotType, ...
            datestr(finneeStc.dataset{ii}.trace{jj}.description.dateOfCreation));
    end
end
fprintf('\n______________________________________________________________________________\n\n')


end

