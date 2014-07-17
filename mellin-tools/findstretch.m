function [ sces, sxcc ] = findstretch( x, s )
%SXCORR  Computes the scale cross-correlation between vectors X and S
%   [SCES SXCC] = FINDSTRETCH(X,S) computes the scale or stretch estimate 
%   between two signals as well as the the scale-invariant correlation 
%   coefficient. The estimates occur in two stages, first estimating the 
%   scale factor through the scale cross-corrleation function. The estimate
%   is then refined in the scale domain, assuming convexity. 
%
%   INPUTS:
%       X: An N-by-K matrix of time-domain signals with N samples,
%          corresponding to K measurements
%       S: An N-by-K matrix of time-domain signals with N samples,
%          corresponding to K measurements
%
%   OUTPUTS:
%    SCES: A K-by-1 vector of scale estimate(s) between X and S
%    SXCC: A K-by-1 vector of scale-invariant correlation coefficients 
%          between X and S
%
%   see also: sxcorr, fmt, fmincon
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


% MEASURE SAMPLE SIZE
N = size(s,1);          % Number of samples
M = nexpsamp(N);        % Number of samples in exponential domain
m = expdoubleaxis(N);   % Axis of exponential domain

% ESTABLISH SCALE-FREQUENCY AXIS
L = log(N)/(M-1);
r = floor(M/2)+1; c = ifftshift((((1:M)-r)/(L*M))).'; 

% FIND SCALE TRANSFORMS AND NORMALIZE
X = fmt(x, 1/2); S = fmt(s, 1/2);
X = X / (norm(X)); S = S / (norm(S));

% MULTIPLY IN THE SCALE DOMAIN
SX = 1/M*bsxfun(@times, conj(X), S);

% FIND THE SCALE CROSS-CORRELATION
sx = ifft(SX); [~, sxi] = max(sx);


% SET OPTIMIZATION OPTIONS 
options = optimset('DiffMaxChange',.01, ...
                   'Algorithm', 'active-set', 'TolX', 1e-12, ...
                   'TolFun', 1e-12, ...
                   'Display', 'notify', 'MaxFunEvals', 1000); 

% IMPROVE ESTIMATE (ASSUMING CONVEXITY AROUND THE MAXIMUM)
sces = fmincon( ...
    @(sces) optfun(X, S, c, sces), ... 
    m(sxi), [], [], [], [], m(sxi)-1/N, m(sxi)+1/N, [], options );

% FIND THE NEW SCALE-INVARIANT CORRELATION COEFFICIENT
sxcc = real(S'*(X.*exp(-1j*2*pi*log(sces)*c)))/norm(X)/norm(S);

end


function [ cx ] = optfun( X, S, c, scale )
    cx = 1-real(S'*(X.*exp(-1j*2*pi*log(scale)*c)));
end