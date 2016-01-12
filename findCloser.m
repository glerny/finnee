function indOut = findCloser(valorIn, vectorIn)
%% DESCRIPTION
%FINDCLOSER is used the closest valor to calorIn in the vector (Xx1 array)
%vectorIn
%   PARAMETERS:
%       valorIn
%           is a number
%
%       vectorIn
%           is a Xx1 array of numbers
%   
%   OUTPUT
%       indOut in the indice in vectorIn where the closest valor to valorIn
%       can be found
%
%% COPYRIGHT
% Copyright 2014-2015 G. Erny (guillaume.erny@finnee.com), FEUP, Porto, Portugal

sTime = tic;
info.type = 'Utility';
info.subFunctionOf = {'getSpectra'};
info.name = 'findCloser';
info.matlabVersion = '8.2.0.701 (R2013b)';
info.version = '20/02/2015_gle01';
info.owner = 'G.Erny @ guillaume.erny@finnee.com';
info.copyright = '';

%% FUNCTION CORE
indMin = find(vectorIn <= valorIn, 1, 'last');
indMax = find(vectorIn >= valorIn, 1, 'first');
if isempty(indMin)
    indOut = indMax;
elseif isempty(indMax)
    indOut = indMin;
elseif abs(vectorIn(indMin) - valorIn) <= abs(vectorIn(indMax) - valorIn)
    indOut = indMin;
else
    indOut = indMax;
end