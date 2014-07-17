function N = nunisamp(M)
%NUNISAMP  Computes the number of samples for an uniformly resampled signal
%   N = NUNISAMP(M) computes the inverse of M = nexpaxis(N), which is the
%   minimum number of samples necessary to exponentially resample a
%   uniformly nyquist sampled signal.
%
%   INPUTS:
%      M: A K-by-1 vector of samples lengths in the exponential domain. 
%
%   OUTPUTS:
%      N: A K-by-1 vector of samples lengths in the uniform domain.
%
%   see also: ifmt, uniaxis, unisamp, nexpsamp
%

% -------------------------------------------------------------------------
% Code written by: Joel B. Harley
% Last updated: July 11, 2014
% -------------------------------------------------------------------------
% If this code is used for a research publication, please cite:
% J.B. Harley, J.M.F. Moura, "Scale transform signal processing for optimal
% ultrasonic temperature compensation," IEEE Transactions on Ultrasonics, 
% Ferroelectrics and Frequency Control, vol. 59, no. 10, October 2012.
% -------------------------------------------------------------------------
%

% CHECK FOR ERRORS
error(nargchk(1, 1, nargin));

% NOTE: This is an approximation but it works extremely well for values of 
%       M less than approx. 90,000,000 samples. 
N = ceil(exp(lambertw(M-1)));

end