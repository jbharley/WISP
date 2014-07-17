function y = sxcorr( x, s, beta, scaleopt, M, itype )
%SXCORR  Computes the scale cross-correlation between vectors X and S
%   Y = SXCORR(X, S, BETA, SCALEOPT, M, ITYPE) computes the scale 
%   cross-correlation between X and S using Mellin real parameter BETA and 
%   interpolation ITYPE.
%
%   INPUTS:
%       X: An N-by-K matrix of time-domain signals with N samples,
%          corresponding to K measurements
%       S: An N-by-K matrix of time-domain signals with N samples,
%          corresponding to K measurements
%    BETA: [OPTIONAL] Real scalar parameter of Mellin transform. By  
%          default, BETA = 1/2.
%SCALEOPT: [OPTIONAL] Type of correlation normalization. Options are
%          'biased', 'unbiased', 'coeff', 'none'. Default is 'none'.
%       M: [OPTIONAL] The scalar number of samples in the Mellin domain.
%          By default, this is approximatly M = N*log(N). See nexpsamp.
%   ITYPE: [OPTIONAL] Type of interpolation. Default is 'spline'.
%          Other choices in 'linear', 'nearest', and 'cubic'
%
%   OUTPUTS:
%       Y: An M-by-K scale cross-correlation signals from X and S. 
%
%   see also: fmt, expsamp, expdoubleaxis
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

% SET DEFAULT MELLIN PARAMETER
if (nargin < 3), beta = 1/2; end;
if (nargin < 4), scaleopt='none'; end;
if (nargin < 5), M = nexpsamp(size(x, 1)); end
if (nargin < 6), itype='spline'; end;

% CHECK IF TRANSPOSE NEEDED
ntrx = 0; if size(x,1) == 1, x = x(:); ntrx = 1; end;
if size(s,1) == 1; s = s(:); end

% CHECK FOR ERRORS
error(nargchk(2, 5, nargin));
if isscalar(M) ~= 1
    error('Error: M should be a scalar value');
end
if size(x,1) ~= size(s,1)
    error('Dimensions of X and S do not match. ');
end

    
    
% NUMBER OF SAMPLES
N = size(x,1);      % number of samples
%M = nexpsamp(N);    % number of exponential domain samples

% EXPONENTIALLY RESAMPLE SIGNALS
xe = expsamp(x, beta, M, itype);
se = expsamp(s, beta, M, itype);

% CORRELATE AND RETURN TO EXPONENTIAL-TIME DOMAIN
Y = bsxfun(@times, conj(fft(xe)), fft(se)); 
if strcmpi(scaleopt, 'biased')
    y = ifft(Y)/M;
elseif strcmpi(scaleopt, 'unbiased')
    y = ifft(Y)/(M-1);
elseif strcmpi(scaleopt, 'coeff')
    y = bsxfun(@times, ifft(Y), 1./(sqrt(sum(abs(xe).^2,1)).*sqrt(sum(abs(se).^2,1))) );
else
    y = ifft(Y);
end
if ntrx; y = y.'; end;


end

