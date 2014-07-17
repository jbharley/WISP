function M = nexpsamp(N)
%NEXPSAMP  Computes number of samples for an exponentially resampled signal
%   M = NEXPSAMP(N) computes the maximum number of samples required to keep
%   the smallest sampling rate in an exponentially resampled signal
%   larger than the original uniform nyquist rate.
%
%   INPUTS:
%      N: A K-by-1 vector of samples lengths in the uniform domain.
%
%   OUTPUTS:
%      M: A K-by-1 vector of samples lengths in the exponential domain. 
%
%   see also: fmt, expaxis, expsamp, nunisamp
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

% COMPUTE THE NUMBER OF NECESSARY SAMPLES
M = ceil(log((N-1)./N.^2) ./ log((N-1)./N));


end