function [ s ] = gaussfilt( N, fc, bw, fs, fnorm )
%GUASSIAN  Compute a Gaussian curve
%   S = GAUSSFILT(N,FC,BW,FS,ENORM) computes an amplitude modulated 
%   Gaussian signal with center frequency FC and bandwith BW
%
%   INPUTS:
%       N: A scalar value denoting the length of the signal
%      FC: A scalar value denoting the center frequency of the signal
%      BW: A scalar value denoting the 3-dB frequency bandwidth of the 
%          Gaussian signal
%      FS: [OPTIONAL] A scalar value denoting the sampling frequency. By 
%          default, FS=1
%   FNORM: [OPTIONAL] If true, Gaussian signal has a maximum amplitude of 
%          one in the frequency domain. If false, Gaussian has a maximum 
%          amplitude of one in the time domain. By default, FNORM = false;
%
%   OUTPUTS:
%       S: The computed Gaussian signal
%

% -------------------------------------------------------------------------
% Code written by: Joel B. Harley
% Last updated: July 16, 2014
% -------------------------------------------------------------------------
%

    % SET DEFAULT PARAMETERS
    if (nargin < 4),    fs = 1; end
    if (nargin < 5), fnorm = false; end
    
    % DEFINE AXES AND PARAMETERS
    t = (0:(N-1)).';
    pw = (1/(2*pi*bw/fs));                 
    mu = N/2;
    sigma = pw * (2*sqrt(2*log(2)));       % Standard deviation
    
    % BUILD GAUSSIAN SIGNAL AND SHIFT TO HAVE ZERO GROUP DELAY
    if fnorm
        s = 2*(1/(2*pi*sigma^2)^(1/2) * exp(-(t - mu).^2./(2*sigma^2))).';
    else
        s = exp(-(t - mu).^2./(2*sigma^2)).';
    end
    s = fshift(ammod(s, fc, fs), N/2).'; 

end

