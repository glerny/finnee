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
    finneeStc.info.parameters.fileID, ...
    finneeStc.info.parameters.folderOut);
% 2. File description
try
fprintf('Original File: %s \n\t(%s  with %s in %s) \n',...
    finneeStc.mzML.fileDescription.sourceFileList{1}.sourceFile{1}.attributes{2}.field{1},...
    finneeStc.mzML.fileDescription.sourceFileList{1}.sourceFile{1}.cvParam{2}.attributes{3}.field{1},...
    finneeStc.mzML.fileDescription.fileContent{1}.cvParam{1}.attributes{3}.field{1},...
    finneeStc.mzML.fileDescription.fileContent{1}.cvParam{2}.attributes{3}.field{1});
catch
end
try
fprintf('Contact Name: \t%s\nContact Organisation: \t%s\n', ...
    finneeStc.mzML.fileDescription.contact{1}.cvParam{1}.attributes{4}.field{1},...
    finneeStc.mzML.fileDescription.contact{1}.cvParam{2}.attributes{4}.field{1});
catch
end
try
fprintf('Date of recording: \t%s\n\n',...
    finneeStc.mzML.run.attributes{1, 3}.field{1});
catch
end
try
fprintf('Sample Name: \t%s\nSample Description: \n', ...
    finneeStc.mzML.sampleList.sample{1}.attributes{2}.field{1})
catch
end
try
fprintf('\t\t%s\n', ...
    finneeStc.mzML.sampleList.sample{1}.userParam{1}.attributes{2}.field{:})
catch
end

% 2. Dataset and traces.
nbrDTS = length(finneeStc.dataset);
fprintf('\nNumber of dataset: \t%d\n', nbrDTS);
for ii = 1:nbrDTS
    fprintf('\t%d %s\n', ii, finneeStc.dataset{ii}.name)
    fprintf('\t\tData format:\t%s\n\t\tDate of creation:\t%s\n',...
        finneeStc.dataset{ii}.description.dataFormat, ...
        datestr(finneeStc.dataset{ii}.dateOfCreation));
    nbrTRC = length(finneeStc.dataset{ii}.trace);
    fprintf('\t\tNumber of traces: \t%d\n\n', nbrTRC);
    
    for jj = 1:nbrTRC
        fprintf('\t\t%d@%d %s\n', jj, ii, ...
            finneeStc.dataset{ii}.trace{jj}.name);
        fprintf('\t\t\tDate of creation:\t%s\n\n',...
            datestr(finneeStc.dataset{ii}.trace{jj}.dateOfCreation));
    end
end
fprintf('\n______________________________________________________________________________\n\n')


end

