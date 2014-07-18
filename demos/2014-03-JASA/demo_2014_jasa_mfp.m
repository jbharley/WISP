% -------------------------------------------------------------------------
% demo_2014_jasa_mfp.m
% -------------------------------------------------------------------------
% This script demonstrates matched field processing (delay-and-sum)
% localization on experimental data. 
%
% Some important variables in the loaded META struct: 
%   META.Rx: Locations of receivers in each trial
%   META.Tx: Locations of transmitters in each trial
%   META.s: Transmitter signal
%   META.x: Reciever signal
%
%
% The script uses one particular functions: FMFP.
% 
% FMFP performs a fast matched field processing approximation on the time 
% data X over frequency samples FN. Note that matched field processing can
% be seen as a generalization of 'delay-and-sum' localization as seen in
% the stuructural health monitoring literautre. The results shown here are
% eseentially from a delay-and-sum algorithm.
%
%
% For more information on the concepts behing this code, see: 
%   Data-driven Matched Field Processing for Lamb Wave 
%   Structural Health Monitoring
%   J.B. Harley, J.M.F. Moura
%   Journal of the Acoustical Society of America
%   vol. 135, no. 3, pp. 1231-1244, March 2014
%

% -------------------------------------------------------------------------
% Copyright (C) 2013-2014  Joel B. Harley
%
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 2 of the License, or (at 
% your option) any later version. You should have received a copy of the 
% GNU General Public License along with this program. If not, see 
% <http://www.gnu.org/licenses/>.
% -------------------------------------------------------------------------
% IF THIS CODE IS USED FOR A RESEARCH PUBLICATION, please cite:
%   J.B. Harley, J.M.F. Moura, "Data-driven matched field processing for 
%   Lamb wave structural health monitoring," Journal of the Acoustical 
%   Society of America, vol. 135, no. 3, March 2014.
% -------------------------------------------------------------------------
% Last updated: July 18, 2014
% -------------------------------------------------------------------------
%

addpath('../../ddmfp-tools')
addpath('../../dsp-tools')

clear;
close all;

% ---------------------------------------------
% USER-INPUTTED INFORMATION
% ---------------------------------------------
% EXPERIMENT INFORMATION
fl = 3;  % 1: one 0.5cm hole 
         % 2: one 0.75cm hole
         % 3: two 0.75cm holes

% SWA INFORMATION
fn = 1:3:120;  % Frequencies 

% LOCALIZATION INFORMATION
W  = 0.75;    % Window size around defect
Nx = 100;     % Number of pixels in x direction
Ny = 100;     % Number of pixels in y direction
Beta = 1/5;   % A factor for improving the speed of FDDMFP -- smaller  
              % number makes FDDMFP faster but more approximate

% WINDOWING INFORMATION
vwin = 2000;  % Velcoity (in m/s) for velocity window

% FILTERING INFORMATION
FC = 300e3;         % Filter center frequency (Hz)
BW = 120e3;         % Filter bandwidth (Hz)
c  = 5.11083e+003;  % Approximate group velocity for 300 kHz


%%
% ---------------------------------------------
% GET DATA
% ---------------------------------------------
load('data_20121029_localize_baseline.mat'); metab = meta; cfgb = cfg;
if fl == 1, load('data_20121029_localize_hole_05cmC.mat'); metas = meta; cfgs = cfg; end
if fl == 2, load('data_20121029_localize_hole_075cmC.mat'); metas = meta; cfgs = cfg; end
if fl == 3, load('data_20121029_localize_2holes_075cmC.mat'); metas = meta; cfgs = cfg; end

%%
% ---------------------------------------------
% GET DISPERSION CURVES
% ---------------------------------------------

% COMPUTE DISTANCES BETWEEN SENSORS
Sx = metas.Sx;  % Hole location
Tx = arrayfun(@(ii) (metab.Tx{ii}.'*ones(1,size(metab.Rx{ii},1))).', 1:numel(metab.Rx), 'UniformOutput', false);  % Transmitter locations
d = diag(dist(cell2mat(metab.Rx), cell2mat(Tx.').'));  % Distance between each reciever and transmitter

% DEFINE NARROW-BAND FILTER 
Q  = size(cell2mat(metab.x.'),1);
sf = gaussfilt( Q, FC, BW, metas.Fs, true );

% PRE-PROCESS SIGNALS
xb = cell2mat(metab.x.');              % Baseline signal
xb(1:110,:) = 0;                       % Remove cross-talk
xb = ifft(bsxfun(@times, fft(sf), fft(xb) ));  % Filter data
xb = velwindow(xb, d, vwin/metab.Fs);  % Velocity-based window
Xb = fft(xb);                          % Fourier transform

xs = cell2mat(metas.x.');              % Signal with target
xs(1:110,:) = 0;                       % Remove cross-talk
xs = ifft(bsxfun(@times, fft(sf), fft(xs) ));  % Filter data
xs = velwindow(xs, d, vwin/metab.Fs);  % Velocity-based window
Xs = fft(xs);                          % Fourier transform

s  = meta.s;                           % Excitation signal
s  = ifft(bsxfun(@times, fft(sf), fft(s, Q) ));  % Filter excitation

%%
% ---------------------------------------------
% PERFORM MATCHED FIELD PROCESSING
% ---------------------------------------------

% BUILD SPATIAL GRID
[ dxy, lx, ly ] = get_grid( metas, cfgs, W, [Nx Ny] );

% PERFORM LOCALIZATION
pxl1 = fmfp( c/metas.Fs, dxy, fn, xb-xs, s, 'beta', Beta, 'model', 'envelope-delay');  % Data-driven matched field processing
pxl1 = pxl1 - min(pxl1); % Normalize -- make minimum value zero
pxl1 = pxl1 / max(pxl1); % Normalize -- make maximum value one
pxl10 = reshape(pxl1(:,1), Nx, Ny);  % Shape into grid


%%
% ---------------------------------------------
% PLOT TIME TRACES AND DISPERSION CURVES
% ---------------------------------------------

% PLOT TIME-DOMAIN DATA
figure(1)
N = metab.N; Fs = metab.Fs; t = 1/Fs:1/Fs:N/Fs;
plot(t, xb)
xlabel('Time')
ylabel('Voltage')
title('Windowed baseline data')

% PLOT LOCALIZATION RESULT
Rx0 = unique(cell2mat(metas.Rx),'rows');
S0 = size(metas.Sx,1);  % Number of scatterers
figure(3)
subplot(1,1+S0,1)
imagesc(lx, ly, 10*log10(pxl10/max(max(pxl10))), [-12 0])
hold on;
plot(Sx(:,1), Sx(:,2), 'kx', 'linewidth', 1.5, 'markersize', 7)
plot(Rx0(:,1), Rx0(:,2), 'ks', 'linewidth', 2, 'markersize', 5)
plot(Sx(:,1), Sx(:,2), 'wx', 'linewidth', 1, 'markersize', 7)
plot(Rx0(:,1), Rx0(:,2), 'ws', 'linewidth', 1, 'markersize', 5)
axis([0 1.22 0 1.22]); axis xy; axis square;
xlabel('Plate width [m]'); ylabel('Plate length [m]');
legend('Hole Center', 'Sensors', 'Location', 'SouthOutside'); legend BOXOFF;
hold off;
title('Delay-and-sum localization')

for s = 1:S0
    subplot(1,1+S0,1+s)
    imagesc(lx, ly, pxl10)
    hold on;
    plot(Sx(:,1) , Sx(:,2) , 'kx', 'linewidth', 2, 'markersize', 10)
    plot(Rx0(:,1), Rx0(:,2), 'ks', 'linewidth', 2, 'markersize', 10)
    plot(Sx(:,1) , Sx(:,2) , 'wx', 'linewidth', 1, 'markersize', 10)
    plot(Rx0(:,1), Rx0(:,2), 'ws', 'linewidth', 1, 'markersize', 10)
    axis([metas.Sx(s,1)-0.03 metas.Sx(s,1)+0.03 metas.Sx(s,2)-0.03 metas.Sx(s,2)+0.03]); axis xy; axis square;
    xlabel('Plate width [m]'); ylabel('Plate length [m]');
    legend('Hole Center', 'Sensors', 'Location', 'SouthOutside'); legend BOXOFF;
    hold off;
    title('Zoomed-in');
end


