function redMS = peakPickingMS_gle1(MS)
%% DESCRIPTION
% 1. INTRODUCTION
% peakPickingMS_gle1 is used to reduced the size of an initial MS scan.
% Local maxima are find and equal to the presence of a MS peak. The average
% position is then determine by fitting each local maxima and it neighbourg
% points with a polynomial of degree 2. 
%
% 2. PARAMETERS:
%   .input: MS is a 2xn array that contain the profile spectrum with in the 
%   first column the m/z increments and in the m/z columns the frequency of 
%   ions detected.
%   .output: redMS id the reduced spectrum with in the first column the
%   accurate mass and in the second the corresponding intensity
%
% 3. EXAMPLES:
%       
% 4. COPYRIGHT
% Copyright 2014-2015 G. Erny (guillaume@fe.up,pt), FEUP, Porto, Portugal

%% CORE OF THE FUNCTION
% Finding non nul local maxima
lmax = [1==2;  (MS(2:end-1,2) >= MS(1:end-2,2) ...
    & MS(2:end-1,2) >= MS(3:end,2)); 1 ==2] & MS(:,2) ~= 0;

% calculate accurate masses with the three most intense
% points
jj = 1;
a = [];
b = [];
c = [];
while 1
    if lmax(jj)
        % We used the analytical solution for the polynomial of degree 2 as
        % with 3 points there is only one solution. This is much faster than
        % a fitting procedure.
        a(end+1) = ((MS(jj-1,2)-MS(jj,2))/(MS(jj-1,1)-MS(jj,1))...
            -(MS(jj-1,2)-MS(jj+1,2))/(MS(jj-1,1)-MS(jj+1,1)))...
            /(MS(jj,1)-MS(jj+1,1));
        b(end+1) = ((MS(jj-1,2)-MS(jj+1,2))-a(end)...
            *(MS(jj-1,1)^2-MS(jj+1,1)^2))...
            /(MS(jj-1,1)-MS(jj+1,1));
        c(end+1) = MS(jj, 2) - a(end)*MS(jj, 1)^2 - b(end)*MS(jj, 1);
        if lmax(jj+1)   % Check if the following data points is also a 
                        % local maxima. If yes, skip it.
            jj = jj + 1;
        end
    end
    jj = jj + 1;
    if jj > length(lmax), break; end
end
redMS = (-b./(2*a))';
redMS(:,2) = round((a'.*redMS(:,1).^2) + (b'.*redMS(:,1)) + c');

% redMS is nan when the peak contain only one point
 redMS(isnan(redMS(:,1)) | isnan(redMS(:,2)), :) = [];


