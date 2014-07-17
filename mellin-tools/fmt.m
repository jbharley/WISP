function X = fmt(x, Beta, M, iType)
%FMT  Compute fast Mellin transform with real parameter BETA
%   XE = FMT(X,BETA,M,ITYPE)  computes the fast Mellin transform of x with
%   real parameter BETA. The Mellin domain is constructed with M samples
%   and is generated using ITYPE interpolation.
%
%   INPUTS:
%       x: An N-by-K matrix of time-domain signals with N samples,
%          corresponding to K measurements
%    BETA: [OPTIONAL] The scalar BETA parameter denotes the amount of gain
%          applied to the signal: t^(BETA) 
%       M: [OPTIONAL] The scalar number of samples in the Mellin domain.
%          By default, this is approximatly M = N*log(N). See nexpsamp.
%   iType: [OPTIONAL] Type of interpolation. Default is 'spline'.
%          Other choices in 'linear', 'nearest', and 'cubic'
%
%   OUTPUTS:
%       X: An M-by-K vector, representing the Mellin transform of x.
%
%   see also: ifmt, expaxis, expsamp, nexpsamp
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


% SET DEFAULT PARAMETERS
if (nargin < 2), Beta = 0; end
if (nargin < 3), M = nexpsamp(size(x, 1)); end
if (nargin < 4), iType = 'spline'; end

% CHECK FOR ERRORS
error(nargchk(1, 4, nargin));
if isscalar(M) ~= 1
    error('Error: M should be a scalar value');
end

% CHECK IF TRANSPOSE NEEDED
ntrx = 0; if size(x,1) == 1; x = x(:); ntrx = 1; end;


% EXPONENTIALLY RESAMPLE SIGNAL
xe = expsamp( x, Beta, M, iType );

% FOURIER TRANSFORM
X = fft(xe);

% TRANSPOSE IF NEEDED
if ntrx; X = X.'; end; 


end

