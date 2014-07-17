function xe = expsamp( x, Beta, M, iType)
%EXPSAMP  Exponentially resamples a signal with attenuation
%   XE = EXPSAMP(X,BETA,M,ITYPE)  exponentially resamples the signal and
%   attentuates by 1/t^(BETA). The exponential-time domain has M samples
%   and is generated using ITYPE interpolation.
%
%   INPUTS:
%       x: An N-by-K matrix of time-domain signals with N samples,
%          corresponding to K measurements
%    BETA: [OPTIONAL] The scalar BETA parameter denotes the amount of gain
%          applied to the signal: t^(BETA) 
%       M: [OPTIONAL] The scalar number of samples in the exponential-time 
%          domain. By default, this is approximatly M = N*log(N). See 
%          nexpsamp.
%   iType: [OPTIONAL] Type of interpolation. Default is 'spline'.
%          Other choices in 'linear', 'nearest', and 'cubic'
%
%   OUTPUTS:
%      xe: Exponentially resampled vector of length M.
%
%   see also: fmt, expaxis, nexpsamp
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


% CHECK IF TRANSPOSE NEEDED
ntrx = 0; if size(x,1) == 1; x = x(:); ntrx = 1; end;

% SET DEFAULT PARAMETERS
if (nargin < 2), Beta = 0; end
if (nargin < 3), M = nexpsamp(size(x, 1)); end
if (nargin < 4), iType = 'spline'; end

% CHECK FOR ERRORS
error(nargchk(1, 4, nargin));
if isscalar(M) ~= 1
    error('Error: M should be a scalar value');
end


% INITIALIZE SIZES
N = size(x, 1);     % Length of each signal
K = size(x, 2);     % Number of signals

% COMPUTE EXPONENTIAL-TIME AXIS
n = (1:N).';
m = expaxis(N, M);

% NOTE: The next line is necessary. If m(end) > n(end) from rounding 
% errors, interp1 will interpret this as an exapolation. By default, the 
% interp1 places a NaN at exapolation points for the 'nearest' and 'linear' 
% settings 
m(end) = n(end);  % Fix rounding error 

% PERFORM INTERPOLATION
xe = interp1(n, x, m, iType);

% ATTENUATE SIGNALS
eBeta = m.^Beta;
xe = xe.*repmat(eBeta,1,K);

% TRANSPOSE IF NEEDED
if ntrx; xe = xe.'; end; 


end

