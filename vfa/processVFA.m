function [llsT1Map] = processVFA(rawVFAdata, TR, flipAngles)
%% PROCESSVFA Calculates T1 maps from VFA data, and AFI B1 data.
% Author: Mathieu Boudreau
%
% Adapted from "main_VFA.m" to turn it into a function.
% Last Modified: November 20th 2013 by Mathieu Boudreau
%
%   ***INPUT***
%   rawVFAdata: 4-D array (X,Y,Z,flipangleindex)
%   TR: Repetition time in seconds
%   flipAngles = 1-D vector with flip angles in degrees. The i'th index
%                of this vector corresponds to the i'th "flipangleindex'
%                of VFAdata.
% 
% Code overview: 
% Note: This code in it's current form uniquely implemements the
% linear-least squared procedure of fitting VFA data.
%
% It also requires B1 maps from the AFI method.

%% VFA Sequence and Image parameters 
% All time sequence parameters are in seconds

acqVFAParam.xDim = length(rawVFAdata(:,1,1));
acqVFAParam.yDim = length(rawVFAdata(1,:,1));
acqVFAParam.numAngles = length(rawVFAdata(1,1,:));

%% T1 calculation from LLS procedure 
% Note that T1 values smaller than 0 s or larger than 4 s are dumped to 
% T1 = 0s.

[llsT1Map, ~, ~, ~] = fitLLSVFAdata(rawVFAdata, flipAngles, acqVFAParam, TR);
llsT1Map(llsT1Map>4)=0;
llsT1Map(llsT1Map<0)=0;
llsT1Map=real(llsT1Map);

end

