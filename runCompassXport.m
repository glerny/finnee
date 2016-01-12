function status = runCompassXport
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

pathCompassXport = ...
    'C:\Program Files (x86)\Bruker Daltonik\CompassXport\CompassXport.exe';
[fileName, pathName] = uigetfile(cd, 'Select the Bruker  file', '*.baf');
targetFile = fullfile(pathName, fileName);
[pathstr,name,~] = fileparts(targetFile);
NewName = fullfile(pathstr, [name, '.mzML']); %dx']);
[fileName, pathName] = uiputfile('*.dx', 'Enter the name of the JCAM file', NewName);
destinationFile = fullfile(pathName, fileName);

command = ['"' pathCompassXport  '"' ' -a ' '"' targetFile '"' ...
    ' -o ' '"' destinationFile '"' ' -mode 2' ' -raw 1'];

status = system(command);
end

