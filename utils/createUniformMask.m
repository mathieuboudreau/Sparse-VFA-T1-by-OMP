function [ mask ] = createUniformMask(centralLines, underSampleFactor, sqDim)
%CREATEUNIFORMMASK Creates a random k-space sampling mask, with a fully acquired
%set of center k-space lines, and uniform probability of undersampling.
% 
% ***Inputs***
% centralLines: Number of central k-space lines
% undersampleFactor: Ratio of total number of lines to number of acquired lines.
% sqDim: Number of voxels in each dimensions (square)
%
% Author: Mathieu Boudreau
% LAST EDITED BY: Mathieu Boudreau
% LAST MODIFIED ON: November 26th 2013

%Pre-allocate array
mask = zeros(1,sqDim);

% Fill center of k-space
mask(1,((sqDim-centralLines)/2+1):(sqDim+centralLines)/2)=1;

% Fill in randomly the rest of k-space until the underSampleFactor is
% achieved
totalLines = sum(mask);
while totalLines<(sqDim/underSampleFactor)
    nextLine = round(1 + (sqDim-1)*rand(1));
    if mask(1,nextLine) ==0
        mask(1,nextLine) = 1;
    end
    totalLines = sum(mask);
end

% Frequency encore lines are fully acquired, copy the acquired phase
% encoded lines for the rest of k-space.
mask = repmat(mask,sqDim,1);

end

