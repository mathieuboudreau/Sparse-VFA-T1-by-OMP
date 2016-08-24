%% OMP-CS VFA T1
% Version: 1.0
%
% DESCRIPTION: This code investigates the effect of Orthogonal Matching  
% Pursuit compressed censing processing applied to noisy 2D VFA data and it's
% effect on T1 values. 
%
% The algorithm uses complete dictionnaries; the dictionnary is NOT sparsified
% using a K-SVD algorithm on training data (but feel free to add this
% feature if you want/need it).
%
% It currently uses a modified Shepp-Logan phantom, but the could should be
% easily modifyable to any T1 maps. The measurements magnitude are 
% calculated from T1 maps, inverse fourier transformed to K-Space and 
% retroactively undersampled prior to any OMP reconstruction.
%
% ***INSTALLATION***
% ************************************************************************* 
% The OMP toolbox (V10) needs to be compiled on your system using MEX. Go
% to ompbox10/readme.txt for more instructions, or go to
% http://www.cs.technion.ac.il/~ronrubin/software.html  to download OMP V10
% *************************************************************************
%
% ***Functions implemented from other open source packages***
% 
% 1. utils/fft2c.m, utils/iff2c.m, utils/genPDF.m and utils/genSampling
% originate from Michael Lustigs Sparse MRI software
% http://www.eecs.berkeley.edu/~mlustig/Software.html
%
% 2. utils/ricernd.m originates from Ged Ridgway's Matlab Rician package
% https://www.mathworks.com/matlabcentral/fileexchange/14237-ricerician-distribution 
%
% 3. ompbox10/ was created by Ron Rubinstein
% http://www.cs.technion.ac.il/~ronrubin/software.html 
%
% Author: Mathieu Boudreau
% LAST EDITED BY: Mathieu Boudreau
% LAST MODIFIED ON: November 26th 2013
%


%% Clear Matlab Session, check dependencies, add paths
%

clear all
close all
clc

addpath([cd '/vfa'])
addpath([cd '/utils'])
addpath([cd '/ompbox10'])

%% Setup Measurements Details
%

param_VFA.TR = 0.015; % in seconds
param_VFA.FA = 2:4:38; % in degrees
param_VFA.sqDim = 256; % Our special Shepp-Logan T1 phantom currently requires 256x256 dimensions

% Noise flag/properties. Noise is added using a Rician distribution.
noiseFlag = 1; % 1 = noisy, 0 = ideal/no noise
SNR = 40;

%% Setup Phantom
%

% Tissue T1 values (seconds)
T1.WM = 0.9;
T1.GM = 1.3;
T1.CSF = 3;

[t1Phantom, maskPhantom]=createModSheppLoganBrainTissueT1Phantom(T1, param_VFA.sqDim); % Includes some simulated lesions!
t1Phantom(:,:,1) = t1Phantom(:,:,1)*100; % M0

imagesc(t1Phantom(:,:,2))
title('Ideal T1 Phantom')
axis image
colorbar

%% Setup Measurements Using T1 map
%

idealdata = calcVFAMeas(t1Phantom,param_VFA.TR,param_VFA.FA);

idealdata(isnan(idealdata)) = 0; % Just in case, but shouldn't happen.

% Add noise?
if noiseFlag == 1
   rawdata = ricernd(idealdata,max(idealdata(param_VFA.sqDim/2,param_VFA.sqDim/2,:))/SNR); % Noise is Ricianly distributed due to magnitude data.
elseif noiseFlag == 0
   rawdata = idealdata;
else
    error('Noise flag not set properly.')    
end

%% Setup k-space mask undersampling properties
%

param_undersampling.undersampleFlag = 1;             % 1 = undersample, 0 = fully sampled
param_undersampling.mask_flag.flag = 'Create';       % Choose between Create or Load
%param_undersampling.mask_flag.file ='filename.mat'; % Set if flag is Load
param_undersampling.underSampleMethod = 'Lustig';    % Choose between Lustig or Uniform; 
                                                     % Both fully samples 16 center lines of k-space. 
                                                     % Lustig is a polynomial density, the polynomial depends on the R factor.
                                                     % Uniform has uniform probability density.
                                                   
Rfactor = [2 4 6]; % K-space acceleration factors. Demo requires three values for the final plot.

%% Launch OMP-Compresse Sensing 
%

% Pre-allocate arrays
llsT1Map=zeros(param_VFA.sqDim,param_VFA.sqDim,length(Rfactor));
recon_data=zeros(param_VFA.sqDim,param_VFA.sqDim,length(param_VFA.FA),length(Rfactor));


for k = 1:length(Rfactor)
    [llsT1Map(:,:,k), recon_data(:,:,:,k)] = processCompSens(rawdata, param_VFA, param_undersampling, Rfactor(k));
end


%% Calculate Normalised Root Mean Square Error of T1 maps, relative to Ideal T1 map.
%

% Pre-allocate arrays
NRMSE_T1 = zeros(1,length(Rfactor));

% Calculate NRMSE of undersampled data
for k = 1:length(Rfactor)
    diff = (llsT1Map(:,:,k).*maskPhantom(:,:) - t1Phantom(:,:,2)).^2/(param_VFA.sqDim^2);
    RMSE_T1 = sqrt(sum(diff(:)));
    NRMSE_T1(k) = RMSE_T1/(max(max(llsT1Map(:,:,k).*maskPhantom(:,:)))-min(min(llsT1Map(:,:,k).*maskPhantom(:,:))));
end

%% Fit VFA T1 of Noisy Fully Acquired Data (R=1)
%

[rawllsT1Map] = processVFA(rawdata(:,:,:,1), param_VFA.TR, param_VFA.FA);

% Calculate NRMSE of fully acquired data
diff = (rawllsT1Map.*maskPhantom(:,:) - t1Phantom(:,:,2)).^2/(param_VFA.sqDim^2);
RMSE_T1 = sqrt(sum(diff(:)));
NRMSE_T1(4) = RMSE_T1/(max(max(rawllsT1Map.*maskPhantom(:,:)))-min(min(rawllsT1Map.*maskPhantom(:,:))));

%% Plot comparison of R=1,2,4,6
%

plotOMPLustigDemo(rawllsT1Map, llsT1Map, maskPhantom, NRMSE_T1)
