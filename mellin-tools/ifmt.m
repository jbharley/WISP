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