function recon_data = ompVFA(numMeas, FA, TR, eps, sqDim, magnImageSpace, phaseImageSpace, rawKspace, mask)
%OMPVFA Reonctructs undersampled VFA data using the OMP algorithm
%   
% Several default values are set within this function, such as the
% plausible T1 range (defaultT1range), the number of atoms in the
% dictionnary (defaultNumAtoms), and the number of non-zero elements in 
% the signal representations.
%
% ***INPUTS*** (The author could certainly shorten this list when he has time)
%
% numMeas: Number of VFA measurements
% FA: VFA excitation flip angles in deg
% TR: Repetition time in seconds
% eps: The convergance cuttoff, variable, depends on images. 
% sqDim: Size of square image
% magnImageSpace: Image of all measurements (magnitude)
% phaseImageSpace: Image of all measurements (phase)
% rawKspace: Complex undersampled k-space 
% mask: K-space undersampling mask.
%
% Author: Mathieu Boudreau
% LAST EDITED BY: Mathieu Boudreau
%

   
%% Create dictionnary
%

defaultT1range = [0.04 4]; %in seconds
defaultNumAtoms = 500;

initialDict = createDictionnary(numMeas,FA,TR, defaultNumAtoms, defaultT1range);

%% Solve OMP-MRI/VFA reconstruction problem
%

% Pre-allocate arrays
resizedMeasurements = zeros(numMeas,sqDim^2);
reconKspace=zeros(sqDim,sqDim,numMeas);
reconMagnImageSpace=zeros(sqDim,sqDim,numMeas);
reconPhaseImageSpace=zeros(sqDim,sqDim,numMeas);

% Start OMP image recon
prevKspace = rawKspace;

for n = 1:1000
        
    count = 1;
    
    % This could certainly be vectorizes
    for k=1:sqDim
       for j=1:sqDim
              for l = 1:numMeas
              resizedMeasurements(l,count) = magnImageSpace(k,j,l);
              end
              count=count+1;
       end
    end

    params.data = resizedMeasurements;
    params.Tdata = 5;               %Number of non-zero elements in the signal representations.
    params.initdict = initialDict;
    
    % Normalize dictionnary signal amplitudes
    params.initdict=params.initdict./repmat(sqrt(sum(params.initdict.^2)),[size(params.initdict,1) 1]);
    
    % Run OMP algorithm
    gamma = omp(params.initdict'*params.data, params.initdict'*params.initdict, params.Tdata);
    
    % Calculate resulting images
    resizedReconImages = params.initdict*full(gamma);
    
    % Reconstruct images into image matrix format
    reconImages = zeros(sqDim,sqDim,numMeas);

    % This could certainly be vectorizes
    count=1;
    for k=1:sqDim
       for j=1:sqDim
              for l = 1:numMeas
              reconImages(k,j,l) = resizedReconImages(l,count);
              end
              count=count+1;
       end
    end
    
    reconImages = abs(reconImages).*exp(1i*phaseImageSpace);
    for k=1:numMeas
            %Fill in original k-space with recon
            reconKspace(:,:,k) = fft2c(reconImages(:,:,k));
            
            reconKspace(:,:,k) = prevKspace(:,:,k).*mask(:,:,k) + reconKspace(:,:,k).*(~mask(:,:,k));
            
            % Transform back to Image domain
            reconMagnImageSpace(:,:,k) = abs(ifft2c(reconKspace(:,:,k)));
            reconPhaseImageSpace(:,:,k) = angle(ifft2c(reconKspace(:,:,k)));
    end
    
    % Calculate error
    calcError = sum(sum(sum( abs(reconMagnImageSpace - magnImageSpace )./(abs (reconMagnImageSpace)) )));
    
    disp(calcError)
    
    % Check if error is below convergence limit set by eps
    if calcError < eps
        recon_data = reconMagnImageSpace;
        break
    else
        magnImageSpace = reconMagnImageSpace;
        phaseImageSpace = reconPhaseImageSpace;
        prevKspace = reconKspace;
    end
    
end

%% If the OMP recon doesn't return recon_data, set as zeros.
%

if ~exist('recon_data','var')
    recon_data = zeros(sqDim,sqDim,numMeas);
end

end

