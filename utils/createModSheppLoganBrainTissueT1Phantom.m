function [t1Phantom, mask] = createModSheppLoganBrainTissueT1Phantom(T1, dimension)
%CREATEMODSHEPPLOGANBRAINTISSUET1PHANTOM Creates a T1 phantom that's a
%little bit more interesting to the author than the conventional 
%Shepp-Logan phantom.
%
% ***INPUTS***
% T1: Structure with the elements WM, GM and CSF containing the T1 values
% in seconds corresponding to these brain tissues.
%
% dimension: Number of element of the (square) image.
%
% Author: Mathieu Boudreau
% LAST EDITED BY: Mathieu Boudreau

WM = T1.WM;
GM = T1.GM;
CSF = T1.CSF;

%Create t1Phantom. "1" in the 3rd dimension is the Magnitude, "2" is
%the T1 value in seconds.

t1Phantom(:,:,1) = phantom(dimension);


% Setup Normal Tissue T1 values
temp_phantom = phantom(dimension)*100;  % Multiplication by 100 needed
% to avoid confusion with T1
% values.

% Set regions to T1
temp_phantom(abs(temp_phantom-20)<0.001) = WM;
temp_phantom((abs(temp_phantom-30)<0.001)|...
    (abs(temp_phantom-40)<0.001)|(abs(temp_phantom-10)<0.001)) = GM;
temp_phantom(abs(temp_phantom-100)<0.001) = CSF;

% Set the center ovals to CSF;
temp_phantom((abs(temp_phantom)<0.0001)&temp_phantom~=0) = WM/2;

t1Phantom(:,:,2) = temp_phantom;
clear temp_phantom


% Create mask
mask = ones(size(t1Phantom,1),size(t1Phantom,2));
mask(t1Phantom(:,:,1)==0) = 0;


% Normalize magnitude of pixels to 1
t1Phantom(:,:,1) = mask;

% Except for the area overlapping zeros, give less tissue there
temp_phantom = mask;
temp_phantom(abs(phantom(dimension)-0.1)<0.001) = 0.5;
% And WM overlapping WM
temp_phantom(abs(phantom(dimension)-0.4)<0.001) = 1.5;

t1Phantom(:,:,1) = temp_phantom;
clear temp_phantom;

% Ensures that no NaN pop up later in the simulation.
t1Phantom(:,:,2) = t1Phantom(:,:,2).*mask;

% Add some fake square lesions
lesionmask = ones(size(t1Phantom,1),size(t1Phantom,2));
lesionmask(72:76, 112:116) = 1.25;
lesionmask(90:94, 140:144) = 1.25;
lesionmask(184:188,84:88) = 1.25;
lesionmask(205:207,135:137) = 1.25;

t1Phantom(:,:,2) = t1Phantom(:,:,2).*lesionmask;
    
end

