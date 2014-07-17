function x = unisamp( xe, Beta, N, iType  )
%UNISAMP  Uniformly resamples (from exponentially sampled) a signal
%   XE = UNISAMP(X,BETA,M,ITYPE)  uniformly resamples a exponentially 
%   sampled signal with attentuation 1/t^(BETA). The uniform-time 
%   domain has N samples and is generated using ITYPE interpolation.
%
%   INPUTS:
%       x: An M-by-K matrix of time-domain signals with M
%          exponential-domian samples, corresponding to K measurements
%    BETA: [OPTIONAL] The scalar BETA parameter denotes the amount of 
%          attenuation applied to each signal: 1 / t^(BETA) 
%       N: [OPTIONAL] The scalar number of samples in the uniform-time 
%          domain. By default, this is approximatly N = nunisamp(M).
%   iType: [OPTIONAL] Type of interpolation. Default is 'spline'.
%          Other choices in 'linear', 'nearest', and 'cubic'
%
%   OUTPUTS:
%      xe: An N-by-K matrix of time-domain signals with N uniform-domian 
%          samples, corresponding to K measurements
%
%   see also: ifmt, uniaxis, nunisamp
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

% CHECK IF TRANSPOSE NEEDED
ntrx = 0; if size(xe,1) == 1; xe = xe(:); ntrx = 1; end;

% SET DEFAULT PARAMETERS
if (nargin < 2), Beta = 0; end
if (nargin < 3), N = nunisamp(size(xe, 1)); end
if (nargin < 4), iType = 'spline'; end

% CHECK FOR ERRORS
error(nargchk(1, 4, nargin));
if isscalar(N) ~= 1
    error('Error: M should be a scalar value');
end



% INITIALIZE SIZES
M = size(xe, 1);     % Length of each signal
K = size(xe, 2);     % Number of signals

% COMPUTE EXPONENTIAL-TIME AXIS
n = (1:N).';
m = expaxis(N, M);

% APPLY GAIN TO SIGNALS
eBeta=m.^(-Beta);
xe=xe.*repmat(eBeta, 1, K);

% PERFORM INTERPOLATION
x = interp1 (m, xe, n, iType);

% TRANSPOSE IF NEEDED
if ntrx; x = x.'; end; 


end

