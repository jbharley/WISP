function y = fshift(x,s)
%FSHIFT  Fractional circular shift
%   Y = FSHIFT(X,S) circularly shifts the elements of vector X by the 
%   possibly non-integer values of s. Shifting is implemented by 
%   multiplying a linear phase term to the Fourier domain of X. 
%
%   INPUTS:
%       X: A Q-by-M matrix of signals with Q samples, corresponding to M
%          measurements
%       S: A scalar value respresnting the factor to shift each signal 
%          by OR a M-by-1 vector of scale factors to shift each signal by
%
%   OUTPUTS:
%       Y: A Q-by-M matrix of shifted signals
% 
%   see also: circshift, fft
%

% -------------------------------------------------------------------------
% Code modified by: Joel B. Harley
% Last updated: July 16, 2014
% -------------------------------------------------------------------------
% This is code is a modified version of the fshift function written by: 
% (c) 2005 Francois Bouffard, All rights reserved.
%          fbouffar@gel.ulaval.ca
%          http://www.mathworks.com/matlabcentral/fileexchange/7886-fshift
%
% From original fshift function: 
% -------------------------------------------------------------------------
% Redistribution and use in source and binary forms, with or without
% modification, are permitted provided that the following conditions are
% met:
% 
%     * Redistributions of source code must retain the above copyright
%       notice, this list of conditions and the following disclaimer.
%     * Redistributions in binary form must reproduce the above copyright
%       notice, this list of conditions and the following disclaimer in
%       the documentation and/or other materials provided with the distribution
% 
% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
% AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
% IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
% ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE
% LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
% CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
% SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
% INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
% CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
% ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
% POSSIBILITY OF SUCH DAMAGE.
%



% CHECK IF TRANSPOSE NEEDED
ntrx = 0; if size(x,1) == 1, x = x(:); ntrx = 1; end;
if size(s,1) == 1; s = s(:); end;

% CHECK FOR ERRORS
if size(x,2) ~= 1 && size(x,2) ~= size(s,1)
    error('Incorrect dimensions for X.');
end

% NUMBER OF SAMPLES
N = size(x,1);   

% DEFINE FREQUENCY VECTOR
r = floor(N/2)+1; f = ifftshift((((1:N)-r)/(N))); 

% DEFINE PHASE VECTOR
p = exp(-1j*2*pi*s*f).';

% PERFORM SHIFT(S)
y = ifft(bsxfun(@times, fft(x), p));
if isreal(x); y = real(y); end;
if ntrx; y = y.'; end; 

end
