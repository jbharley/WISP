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
% Code written by: Joel B. Harley
% Last updated: July 11, 2014
% -------------------------------------------------------------------------
% If this code is used for a research publication, please cite:
% J.B. Harley, J.M.F. Moura, "Scale transform signal processing for optimal
% ultrasonic temperature compensation," IEEE Transactions on Ultrasonics, 
% Ferroelectrics and Frequency Control, vol. 59, no. 10, October 2012.
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

