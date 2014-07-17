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
% Copyright (C) 2012,2014  Joel B. Harley
%
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 2 of the License, or (at 
% your option) any later version. You should have received a copy of the 
% GNU General Public License along with this program. If not, see 
% <http://www.gnu.org/licenses/>.
% ------------------------------------------------------------------------- 
% IF THIS CODE IS USED FOR A RESEARCH PUBLICATION, please cite:
%   J.B. Harley, J.M.F. Moura, "Sparse Recovery of the Multimodal and 
%   Dispersive Characteristics of Lamb Waves," Journal of the Acoustical
%   Society of America, vol. 133, no. 5, May 2013.
% -------------------------------------------------------------------------
% Last updated: July 16, 2014
% -------------------------------------------------------------------------
%

% CHECK FOR ERRORS
error(nargchk(1, 1, nargin));

% NOTE: This is an approximation but it works extremely well for values of 
%       M less than approx. 90,000,000 samples. 
N = ceil(exp(lambertw(M-1)));

end