% -------------------------------------------------------------------------
% demo_2014_jasa.m
% -------------------------------------------------------------------------
% This script demonstrates data-driven matched field processing
% localization on experimental data. 
%
% Some important variables in the loaded META struct: 
%   META.Rx: Locations of receivers in each trial
%   META.Tx: Locations of transmitters in each trial
%   META.s: Transmitter signal
%   META.x: Reciever signal
%
%
% The script uses two particular functions: SWA and FDDMFP.
%
% SWA performs sparse wavenumber analysis on the time data X over frequency
% samples FN and over a wavenumbers K. Sparse wavenumber analysis
% can be performed using one of two methods currently: basis pursuit 
% denoising ('bp') or orthogonal matching pursuit* ('omp'). 
%
% * Note that I am currently preparing a JASA Express letter discussing the 
% use of orthogonal matching pursuit (omp) for sparse wavenumber analysis  
% and its tradeoffs with basis pursuit denoising. Paper title can be found 
% at the end of these comments. 
% 
% FDDMFP performs a fast data-driven matched field processing approximation
% on the time data X % over frequency samples FN with the aid of dispersion 
% curves V obtained from SWA. 
%
%
% For more information on the concepts behing this code, see: 
%   Sparse Recovery of the Multimodal and Dispersive  
%   Characteristics of Lamb Waves
%   J.B. Harley, J.M.F. Moura
%   Journal of the Acoustical Society of America
%   vol. 133, no. 5, pp. 2732-2745, May 2013
%
%   Data-driven Matched Field Processing for Lamb Wave 
%   Structural Health Monitoring
%   J.B. Harley, J.M.F. Moura
%   Journal of the Acoustical Society of America
%   vol. 135, no. 3, pp. 1231-1244, March 2014
%
% The analysis and comparison of a orthogonal matching pursuit 
% implementation of sparse wavenumber analysis (used this code) and its 
% a basis pursuit denoising implementation (used in above papers) will 
% soon be submitted in: 
%   Fast Dispersion Curve Recovery with Orthogonal Matching Pursuit
%   J.B. Harley, J.M.F. Moura
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
% Last updated: July 16, 2014
% -------------------------------------------------------------------------
%

addpath('../../ddmfp-tools')
addpath('../../dsp-tools')

clear;
% ---------------------------------------------
% USER-INPUTTED INFORMATION
% ---------------------------------------------
% EXPERIMENT INFORMATION
fl = 3;  % 1: one 0.5cm hole 
         % 2: one 0.75cm hole
         % 3: two 0.75cm holes

% SWA INFORMATION
fn = 100:100:800;   % Frequencies 
k  = (2:2:1000).';  % Wavenumbers for dispersion curves

% LOCALIZATION INFORMATION
W  = 1;       % Window size around defect
Nx = 500;     % Number of pixels in x direction
Ny = 500;     % Number of pixels in y direction
Beta = 1/10;  % A factor for improving the speed of FDDMFP -- smaller  
              % number makes FDDMFP faster but more approximate

% WINDOWING INFORMATION
vwin = 2000;  % Velcoity (in m/s) for velocity window
              
% MATCHED FIELD PROCESSORS
incoherentProc = false;  % Incoherent matched field processor

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

% PRE-PROCESS SIGNALS
xb = cell2mat(metab.x.');              % Baseline signal
xb(1:110,:) = 0;                       % Remove cross-talk
xb = velwindow(xb, d, vwin/metab.Fs);  % Velocity-based window
Xb = fft(xb);                          % Fourier transform

xs = cell2mat(metas.x.');              % Signal with target
xs(1:110,:) = 0;                       % Remove cross-talk
xs = velwindow(xs, d, vwin/metab.Fs);  % Velocity-based window
Xs = fft(xs);                          % Fourier transform

% GET DISPERSION CURVES
V = sparse(numel(k), size(Xb,1));
V(:,fn) = swa( k, d, Xb(fn,:), 2, 'method', 'omp' );


%%
% ---------------------------------------------
% PERFORM MATCHED FIELD PROCESSING
% ---------------------------------------------

% BUILD SPATIAL GRID
[ dxy, lx, ly ] = get_grid( metas, cfgs, W, [Nx Ny] );

% PERFORM LOCALIZATION
pxl1 = fddmfp( k, dxy, fn, xb-xs, V, 'beta', Beta, 'incoherent', incoherentProc );  % Data-driven matched field processing
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

% PLOT RECOVERED DISPERSION CURVES
figure(2)
f = linspace(0, Fs, N);
imagesc(f(fn)/1000, k, 10*log10(abs(V(:,fn))/max(max(abs(V(:,fn))))), [-25 0])
axis xy;
ylabel('Wavenumber')
xlabel('Frequency [kHz]')
a = colorbar;
xlabel(a, '[dB]')
title('Recovered dispersion curve (at select frequencies)')

% PLOT LOCALIZATION RESULT
Rx0 = unique(cell2mat(metas.Rx),'rows');
S0 = size(metas.Sx,1);  % Number of scatterers
figure(3)
subplot(1,1+S0,1)
imagesc(lx, ly, pxl10)
hold on;
plot(Sx(:,1), Sx(:,2), 'kx', 'linewidth', 1.5, 'markersize', 7)
plot(Rx0(:,1), Rx0(:,2), 'ks', 'linewidth', 2, 'markersize', 5)
plot(Sx(:,1), Sx(:,2), 'wx', 'linewidth', 1, 'markersize', 7)
plot(Rx0(:,1), Rx0(:,2), 'ws', 'linewidth', 1, 'markersize', 5)
axis([0 1.22 0 1.22]); axis xy; axis square;
xlabel('Plate width [m]'); ylabel('Plate length [m]');
legend('Hole Center', 'Sensors', 'Location', 'SouthOutside'); legend BOXOFF;
hold off;
title('Data-driven matched field processing')

for s = 1:S0
    subplot(1,1+S0,1+s)
    imagesc(lx, ly, pxl10)
    hold on;
    plot(Sx(:,1), Sx(:,2), 'kx', 'linewidth', 2, 'markersize', 10)
    plot(Rx0(:,1), Rx0(:,2), 'ks', 'linewidth', 2, 'markersize', 10)
    plot(Sx(:,1), Sx(:,2), 'wx', 'linewidth', 1, 'markersize', 10)
    plot(Rx0(:,1), Rx0(:,2), 'ws', 'linewidth', 1, 'markersize', 10)
    axis([metas.Sx(s,1)-0.03 metas.Sx(s,1)+0.03 metas.Sx(s,2)-0.03 metas.Sx(s,2)+0.03]); axis xy; axis square;
    xlabel('Plate width [m]'); ylabel('Plate length [m]');
    legend('Hole Center', 'Sensors', 'Location', 'SouthOutside'); legend BOXOFF;
    hold off;
    title('Zoomed-in');
end


