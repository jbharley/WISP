function m = expdoubleaxis (N, M)
%EXPAXIS  Constructs an scale-symmetric exponential-time domain axis
%   m = EXPAXIS(N,M)  constructs an M samples exponential-time axis from a
%   N samples uniform-time axis. This axis is specfically for scale
%   cross-correlated signals. The axis is expontially circularly symmetric 
%   around sample 1.
%
%   INPUTS:
%       N: The scalar number of samples in a uniformly sampled signal
%       M: [OPTIONAL] The scalar number of samples in the desired 
%          exponentially sampled signal. By default, this is approximatly 
%          M = N*log(N). See nexpsamp.
%
%   OUTPUTS:
%       m: An M-by-1 vector representing an exponentially sampled 
%          cross-correlation axis
%
%   see also: expaxis, fmt, expsamp, nexpsamp, sxcorr
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
%   J.B. Harley, J.M.F. Moura, "Scale transform signal processing for 
%   optimal ultrasonic temperature compensation," IEEE Transactions on  
%   Ultrasonics, Ferroelectrics and Frequency Control, vol. 59, no. 10, 
%   October 2012.
% -------------------------------------------------------------------------
% Last updated: July 16, 2014
% -------------------------------------------------------------------------
%


% SET DEFAULT PARAMETERS
if (nargin < 2)
    % Maximum number of samples required to keep the smallest sampling
    % rate larger than the original uniform nyquist rate
    M = nexpsamp(N);
end

% CHECK FOR ERRORS
error(nargchk(1, 2, nargin));
if ~isscalar(N) || ~isscalar(M)
    error('Error: N and M should be scalars.');
end


% COMPUTE SYMMETRIC EXPONENTIAL-TIME AXIS
L = log(N)/(M-1);
r = floor(M/2)+1; c = ifftshift((1:M)-r)*L; 

% CONVERT TO SCALE AXIS
m = exp(c).';


end




