function [VFASignal] = calcVFAMeas(t1Phantom,TR,FA)
%CALCVFAMEAS Calculates 2D VFA measurement maps from a T1 map.
% 
% ***Inputs***
% t1Phantom: 2D T1 map in seconds
% TR: VFA repetition time in seconds
% FA: VFA flip angles in degrees.
%
% Author: Mathieu Boudreau
% LAST EDITED BY: Mathieu Boudreau
    
    VFASignal = zeros(size(t1Phantom,1),size(t1Phantom,2), length(FA));

    for k = 1:length(FA)
        VFASignal(:,:,k) = t1Phantom(:,:,1).*((1-exp(-TR./t1Phantom(:,:,2)))).*sind(FA(k))./(1-exp(-TR./t1Phantom(:,:,2)).*cosd(FA(k))); 
    end


end

