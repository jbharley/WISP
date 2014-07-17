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

% COMPUTE THE NUMBER OF NECESSARY SAMPLES
M = ceil(log((N-1)./N.^2) ./ log((N-1)./N));


end