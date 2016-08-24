function [llsT1Map, llsEquiliMagnMap, xVFAFitTerm, yVFAFitTerm] = fitLLSVFAdata(rawVFAdata, flipAngle, acqVFAParam, TR)
%FITLLSVFA Fits T1 from 2D VFA data using the linear y = mx+b form of the 
%VFA equation.
% 
% This function output the T1 maps fitted from the data using a linear
% least squared algorithm.
%
% Note: A bit slow due to nested for-loops, someone should fix that 
% someday... 
%
% Author: Mathieu Boudreau
% LAST EDITED BY: Mathieu Boudreau

%% Reorder VFA equation to y=mx+b form.
%

% Pre-allocate array
xVFAFitTerm = zeros(acqVFAParam.xDim, acqVFAParam.yDim, 2);
yVFAFitTerm = zeros(acqVFAParam.xDim, acqVFAParam.yDim, 2);

VFAdata = double(rawVFAdata);
for k=1:length(flipAngle)
    xVFAFitTerm(:,:,k) = VFAdata(:,:,k)/tand(flipAngle(k)); % x = Signal/tan(alpha)
    yVFAFitTerm(:,:,k) = VFAdata(:,:,k)/sind(flipAngle(k)); % y = Signal/sin(alpha)
end

%% Fit for y = mx + b
% Note that: fitParamVFA(:,1) = m; fitParamVFA(:,2) = b

% Pre-allocate array
fitParamVFA = zeros(acqVFAParam.xDim, acqVFAParam.yDim, 2);

for x = 1:acqVFAParam.xDim
    disp(x)
    for y = 1:acqVFAParam.yDim
        fitParamVFA(x,y,1:2) = polyfit(squeeze(xVFAFitTerm(x,y,:)), squeeze(yVFAFitTerm(x,y,:)), 1);
    end
end

%% Extract T1 from fit
%

% m = exp(-TR/T1), solve for T1
llsT1Map = -TR./log(fitParamVFA(:,:,1));

% b = M0(1-exp(-TR/T1)), solve for M0
llsEquiliMagnMap = fitParamVFA(:,:,2)./(1-exp(TR./llsT1Map));

end

