function [llsT1Map, recon_data] = processCompSens(rawdata, param_VFA, param_undersampling, Rfactor)
%PROCESSCOMPSENS Retroactively undersamples k-space of VFA data, then
%reconstructs measurement images using OMP.
%
% The variable "eps" in the section "Reconstruct image using OMP" is the
% convergence cutoff criteria for the OMP reconstruction. Any lower than 
% the default for the mainDemo will take a long time through and won't give
% that great improvements. 
%
% ***INPUT***
% rawdata: Array of raw VFA measurement images.
% param_VFA: Structure containing VFA measurement details.
% param_undersampling: Structure containing undersampling details.
% Rfactor: Positive value that defines the acceleration factor
%
% Author: Mathieu Boudreau
% LAST EDITED BY: Mathieu Boudreau
% LAST MODIFIED ON: November 26th 2013
%

%% Deconstruct structures
%

% param_VFA
TR = param_VFA.TR;
FA = param_VFA.FA;
sqDim = param_VFA.sqDim;
numMeas = length(FA);


% param_undersampling
undersampleFlag = param_undersampling.undersampleFlag;
mask_flag = param_undersampling.mask_flag;
underSampleMethod = param_undersampling.underSampleMethod;

%% Setup k-space mask
%

% Pre-allocate arrays
rawKspace = zeros(size(rawdata,1),size(rawdata,2), size(rawdata,3));
magnImageSpace = zeros(size(rawdata,1),size(rawdata,2), size(rawdata,3));
phaseImageSpace = zeros(size(rawdata,1),size(rawdata,2), size(rawdata,3));
mask=zeros(size(rawdata,1),size(rawdata,2), size(rawdata,3));

if undersampleFlag ==0
        mask = ones(sqDim,sqDim,numMeas);

elseif undersampleFlag == 1
    switch mask_flag.flag
        case 'Load'
            load(mask_flag.file)
        case 'Create'
            switch underSampleMethod                
                case 'Lustig'
                    % Pre-allocate arrays
                    pdfTemp = zeros(size(rawdata,1),size(rawdata,3));
                    maskTemp = zeros(size(rawdata,1),size(rawdata,3));
                    pdf = zeros(size(rawdata,1),size(rawdata,2),size(rawdata,3));

                    if abs(Rfactor-2)<0.001
                        polyNom = 2;
                    elseif abs(Rfactor-3)<0.001
                        polyNom = 3;
                    elseif abs(Rfactor-4)<0.001
                        polyNom = 4;
                    else
                        polyNom = 12;
                    end
                    for k = 1:numMeas
                        pdfTemp(:,k)= genPDF([1,sqDim],polyNom,1/Rfactor,2,0.16,0);
                        maskTemp(:,k) = genSampling(pdfTemp(:,k),100,1);

                        pdf(:,:,k) = repmat(pdfTemp(:,k)',sqDim,1);
                        mask(:,:,k) = repmat(maskTemp(:,k)',sqDim,1);
                    end
                case 'Uniform'
                  rng('shuffle');
                  centralLines = 16;                  
                  for k = 1:numMeas
                      mask(:,:,k) = createUniformMask(centralLines, Rfactor, sqDim);
                  end
            end
    end
    
    disp('The r factor is: ')
    disp((sqDim^2)/(sum(mask(:))/numMeas))
    for k=1:numMeas
            %Undersample k-space
            rawKspace(:,:,k) = fft2c(rawdata(:,:,k)).*mask(:,:,k);
            % Transform back to Image domain
            magnImageSpace(:,:,k) = abs(ifft2c(rawKspace(:,:,k)));
            phaseImageSpace(:,:,k) = angle(ifft2c(rawKspace(:,:,k)));
    end
end

%% Reconstruct image using OMP
%
eps = 2000; % *** Convergence criteria, lower means longer processing, this is a value that the author found was decent through purely manual labour ***

recon_data = ompVFA(numMeas, FA, TR, eps, sqDim, magnImageSpace, phaseImageSpace, rawKspace, mask);

%% Process VFA
%

[llsT1Map] = processVFA(recon_data, TR, FA);
end

