function [ d, lx, ly ] = get_grid( meta, cfg, dW, N )
%GET_GRID  Build localization from experimental meta data
%   [D LX LY] = GET_GRID( META, CFG, DW, N ) uses meta data from
%   experiments to generate a grid of distances associated with each
%   recieved measurement
%
%   INPUTS:   
%    META: A struct of meta data from the extra_data function
%     CFG: A struct of configuration data from the extra_data function
%      DW: A scalar value for the minimum spatial window around a target. 
%          The window never exceeds the boundaries of the experiment. 
%          If there are multiple targets, the window center is their 
%          mean position. If there are no targets, the window center is 
%          the medium's center. 
%       N: A scalar or 2-by-1 vector [Nx Ny] representing the number of 
%          grid points on the x- and y-axis. 
%
%   OUTPUTS:
%       D: An Nx*Ny-by-M matrix of distances corresponding to each
%          recieved measurement
%      LX: An Nx-by-1 vector representing the x-axis of the grid
%      LY: An Ny-by-1 vector representing the y-axis of the grid
%
%   NOTE:
%     This function currently does not work for multiple simultaneous
%     sources
%

% -------------------------------------------------------------------------
% Copyright (C) 2014  Joel B. Harley
%
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 2 of the License, or (at 
% your option) any later version. You should have received a copy of the 
% GNU General Public License along with this program. If not, see 
% <http://www.gnu.org/licenses/>.
% -------------------------------------------------------------------------
% Last updated: July 16, 2014
% -------------------------------------------------------------------------
%

    % CHECK INPUTS
    if numel(N) > 2, error('numel(N) should be less than or equal to 2.'); end
    if numel(N) == 1, N = [N N]; end

    % GET BOUNDS OF EXPERIMENT
    XB = [ min(min(cfg.bound(:,1:2))) max(max(cfg.bound(:,1:2))) ];  % X-bounds
    YB = [ min(min(cfg.bound(:,3:4))) max(max(cfg.bound(:,3:4))) ];  % Y-bounds

    % CHOOSE GRID CENTER BASED ON TARGET LOCATION
    if iscell(meta.Sx), Sx = unique(cell2mat(meta.Sx),'rows'); else Sx = meta.Sx;  end
    if isempty(Sx), Sx = [mean(XB) mean(YB)]; else Sx = mean(Sx,1); end
    
    % GET BOUNDS OF GRID
    Xmin = max(XB(1), Sx(1) - dW/2); Xmax = min(XB(2), Sx(1) + dW/2);  % Window X-bounds
    Ymin = max(XB(1), Sx(2) - dW/2); Ymax = min(YB(2), Sx(2) + dW/2);  % Window Y-bounds

    % BUILD GRID
    lx  = linspace(Xmin, Xmax, N(1));  % X-axis
    ly  = linspace(Ymin, Ymax, N(2));  % Y-axis     
    [Lx Ly] = meshgrid(lx,ly);         % Separate axes
    Gxy = [reshape(Lx, N(1)*N(2), 1) reshape(Ly, N(1)*N(2), 1)];  % Define grid
    
    % DEFINE DISTANCES 
    Rx = cell2mat(meta.Rx);  % Reciever locations
    Tx = cell2mat(arrayfun(@(ii) (meta.Tx{ii}.'*ones(1,size(meta.Rx{ii},1))).', 1:numel(meta.Rx), 'UniformOutput', false).');  % Transmitter locations
    if isempty(Tx)
        d = dist(Rx, Gxy.');                    % Active source
    else
        d = dist(Rx, Gxy.') + dist(Tx, Gxy.');  % Passive source
    end
    
    
end

