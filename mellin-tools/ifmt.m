function x = ifmt(X, Beta, N, iType)
%IFMT  Compute inverse fast Mellin transform with real parameter BETA
%   XE = IFMT(X,BETA,M,ITYPE)  computes the fast Mellin transform of x with
%   real parameter BETA. The Mellin domain is constructed with M samples
%   and is generated using ITYPE interpolation.
%
%   INPUTS:
%       X: An M-by-K matrix of time-domain signals with M samples,
%          corresponding to K measurements
%    BETA: [OPTIONAL] The scalar BETA parameter denotes the amount of 
%          attenuation applied to the signal: 1 / t^(BETA) 
%       N: [OPTIONAL] The scalar number of samples in the uniform-time
%          domain. By default, this is approximatly M = nunisamp(M).
%   iType: [OPTIONAL] Type of interpolation. Default is 'spline'.
%          Other choices in 'linear', 'nearest', and 'cubic'
%
%   OUTPUTS:
%       x: An N-by-K vector, representing the inverse Mellin transform 
%          of X.
%
%   see also: fmt, expsamp, unisamp, nunisamp
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
if (nargin < 3), N = nunisamp(size(X, 1)); end
if (nargin < 4), iType = 'spline'; end

% CHECK FOR ERRORS
error(nargchk(1, 4, nargin));
if isscalar(N) ~= 1
    error('Error: N should be a scalar value');
end

% CHECK IF TRANSPOSE NEEDED
ntrx = 0; if size(X,1) == 1; X = X(:); ntrx = 1; end;


% INVERSE FOURIER TRANSFORM
xe = ifft(X);

% RETURN TO A UNIFORMLY SAMPLED SIGNAL
x = unisamp( xe, Beta, N, iType );

% TRANSPOSE IF NEEDED
if ntrx; x = x.'; end; 


end