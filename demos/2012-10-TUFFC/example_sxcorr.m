% -------------------------------------------------------------------------
% example_sxcorr.m
% -------------------------------------------------------------------------
% In this script, we show examples of two functions used to estimate
% the scale factor between two signals and estimate their scale 
% invariant correlation coefficient. The two functions are SXCORR and
% FINDSTRETCH.
%
% SXCORR computes the discrete-time scale cross-correlation, the maximum of
% which is the scale-invariant correlation coefficient and the argument of
% the maximum is the scale estimate. The problem with SXCORR is that its
% resolution is limited by discrete nature of the signal. SXCORR has a
% resolution of approximately 1/N, where N is the length of the signals of
% interest. This approach is known as Scale-Invariant Correlation (SIC) in
% my paper (see citation below)
%
% FINDSTRETCH fixes this issue by using a iterative opimization approach to
% find the exact scale factor. SXCORR is first used to compute the inital
% estimate of the scale factor. An interior point algorithm is then used to
% refine that estimate in the scale domain. This algorithm is more accurate
% at the cost of algorithmic complexity. This approach is known as 
% Scale-Invariant Correlation / Iterative Scale Transform (SIC/IST) in my
% paper (see citation below)
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

addpath('../../mellin-tools/')
addpath('../../dsp-tools/')


clear;

% DEFINE SOME VARIABLES
N = 1000;               % Length of timm axis

dly = 200;              % Delay of Gaussian pulse
fc = .2;                % Center frequency of Gaussian pulse
bw = .2;                % Bandwidth of Gaussian pulse
stretch = 2.12345;      % Specify stretch factor

fprintf('\n')
fprintf('-----------------------------------------------------------------\n')
fprintf(' This is an example of SXCORR and FINDSTRETCH functions \n');
fprintf('-----------------------------------------------------------------\n')
fprintf('True stretch factor = %f \n\n', stretch)


% ----------------------------------------------------------------------

% Time axis
n = 1:N;        

% Axis for the scale-cross correlation result
m = expdoubleaxis(N);


% Generate gaussian pulse using builtin gauspuls equation...
s = fshift(gaussfilt( N, fc, .025, 1, false), dly);

% Stretch signal using circscale function (my function) 
x = circscale(s, stretch);

% Compute scale cross-correlation between s and x for plotting
sx = sxcorr(s,x, 1/2, 'coeff');

% PLOT
figure(1)
subplot(311)
plot(n,s)
xlabel('Time'); ylabel('Amplitude');
title('Baseline signal')
subplot(312)
plot(n, x)
xlabel('Time'); ylabel('Amplitude');
title('Test signal')
subplot(313)
plot(fftshift((m)), fftshift(sx))
xlabel('Scale factor'); title('Scale cross-correlation');



% ----------------------------------------------------------------------
% USING SXCORR
% ----------------------------------------------------------------------

% Compute SXCORR (Beta = 1/2 specifies the scale transform) 
%   ('coeff' specifies a normalization with a maximum of 1)
sx = sxcorr(s,x, 1/2, 'coeff');

% Maximum value is scale-invariant correlation coeffient 
% Index of maximum value provides scale estimate
[sxx sxi] = max(sx);

% Note that the scale-invariant correlation coefficient should be 1 but is
% not due to resolution constraints (since the signals are discrete)

fprintf('SXCORR      : Scale estimate = %f\n', m(sxi));
fprintf('SXCORR      : Scale-invariant correlation coeff. = %f\n', sxx);

% Scale s back based on the estimated scale value
s0 = circscale(x, 1./m(sxi));

% PLOT
figure(2)
subplot(211)
plot(n, s, 'linewidth', 2); hold on;
plot(n, s0, '--r', 'linewidth', 2); hold off;
xlabel('Time'); ylabel('Amplitude');
title('SXCORR results')
legend('Baseline signal', 'Scaled test signal')
subplot(212)
plot(s0-s)
xlabel('Time'); ylabel('Signal difference (error)');


% ----------------------------------------------------------------------
% USING FINDSTRETCH
% ----------------------------------------------------------------------

% Using the FINDSTETCH function
[ sces sxcc ] = findstretch( s, x );

fprintf('FINDSTRETCH : The scale estimate = %f\n', sces);
fprintf('FINDSTRETCH : Scale-invariant correlation coeff. = %f\n', sxcc);

% Scale s back based on the estimated scale value (from FINDSTETCH)
s1 = circscale(x, 1./sces);

% PLOT
figure(3)
subplot(211)
plot(n, s, 'linewidth', 2); hold on;
plot(n, s1, '--r', 'linewidth', 2); hold off;
xlabel('Time'); ylabel('Amplitude');
title('FINDSTRETCH results')
legend('Baseline signal', 'Scaled test signal')
subplot(212)
plot(s1-s)
xlabel('Time'); ylabel('Signal difference (error)');


