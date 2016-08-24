function [dictionnary] = createDictionnary(numMeas, FA, TR, numAtoms, T1range)
%CREATEDICTIONNARY Creates a dictionnary of VFA signal atoms 
%
% ***INPUTS***
% numMeas: Number of VFA measurements
% FA: VFA excitation flip angles in degrees
% TR: Repetition times in seconds
% numAtoms: Total number of atoms
% T1range: Range of T1 values of the dictionnary
%
% Author: Mathieu Boudreau
% LAST EDITED BY: Mathieu Boudreau
%

% Pre-allocate arrays
dictionnary = zeros(numMeas,numAtoms);

T1 = linspace(T1range(1),T1range(2),numAtoms);

% Calculate dictionnary
for k=1:numMeas
    dictionnary(k,:) = 1*((1-exp(-TR./T1)).*sind(FA(k)))./(1-exp(-TR./T1).*cosd(FA(k))); % VFA steady state signal equation
end

end

