function y = circscale( x, s, beta )
%CIRCSCALE Circularly scales a signal using the Mellin transform
%   Y = CIRCSCALE(X,S,BETA) circularly scales the elements of vector X by  
%   the non-integer values of s. Scaling is implemented by multiplying a
%   linear phase term to the Mellin domain of X. BETA determines the real 
%   part of the Mellin transform and the attentuation on X. 
%
%   INPUTS:
%       X: An N-by-K matrix of time-domain signals with N samples,
%          corresponding to K measurements
%       S: A scalar value respresnting the factor to scale, or stretch, 
%          each signal by OR a K-by-1 vector of scale factors to scale 
%          each signal by
%    BETA: [OPTIONAL] Real scalar parameter of Mellin transform. By  
%          default, BETA = 0.
%
%   OUTPUTS:
%       Y: An N-by-K matrix of scaled, or stretched,  time-domain signals 
%          with N samples, corresponding to K measurements
%
%   see also: fmt, expsamp, circshift
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
if (nargin < 3), beta=0; end

% CHECK IF TRANSPOSE NEEDED
ntrx = 0; if size(x,1) == 1; x = x(:); ntrx = 1; end;
if size(s,1) == 1; s = s(:); end;

% CHECK FOR ERRORS
if size(x,2) ~= 1 && size(x,2) ~= size(s,1)
    error('Incorrect dimensions for X.');
end



% NUMBER OF SAMPLES
N = size(x,1);      % number of samples
M = nexpsamp(N);    % number of exponential domain samples

% ESTABLISH SCALE-FREQUENCY AXIS
L = log(N)/(M-1);
r = floor(M/2)+1; c = ifftshift((((1:M)-r)/(L*M))); 

% DEFINE PHASE VECTOR
p = exp(-1j*2*pi*log(s)*c).';

% PERFORM SCALE
y = ifmt(bsxfun(@times, fmt(x, beta), p), beta);
if isreal(x); y = real(y); end;
if ntrx; y = y.'; end; 


end

