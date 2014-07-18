function X = sws( k, d, V )
%SWS  Sparse wavenumber synthesis
%   [VS Y A Dr] = SWA( K, D, V ) synthesizes raw wave data from its 
%   sparse frequency-wavenumber representation
%
%   INPUTS:       
%       K: An N-by-1 vector of wavenumbers
%       D: An M-by-1 vector of distances
%       V: An N-by-Q matrix representing the wave's sparse 
%          frequency-wavenumber representation
%
%   OUTPUTS:
%       X: A Q-by-M matrix of synthesized data in the frequency domain
%
%   see also: swa, ddmfp, fddmfp
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
%   J.B. Harley, J.M.F. Moura, "Sparse Recovery of the Multimodal and 
%   Dispersive Characteristics of Lamb Waves," Journal of the Acoustical
%   Society of America, vol. 133, no. 5, May 2013.
% -------------------------------------------------------------------------
% Last updated: July 16, 2014
% -------------------------------------------------------------------------
%

    % CHECK INPUTS 
    error(nargchk(3, 3, nargin));
    
    % FORCE COLUMN VECTORS
    k = k(:);  % Wavenumber vector
    d = d(:);  % Distance vector
    
    % DEFINE LENGTHS
    M = size(d, 1);  % Number of signals to synthesize
    N = size(k, 1);  % Number of wavenumbers
    
    % COMPUTE PROPOGATION FRAME
    A  = exp(-1j*d*k.');                % Complex exponentials
    Dr = spdiags(1./sqrt(d), 0, M, M);  % Distance attenuation
    Dk = spdiags(1./sqrt(k), 0, N, N);  % Wavenumber attenutation 
    Phi = Dr*A*Dk;                      % Progagation frame 
                                        
    % SYNTHESIZE SIGNALS
    X = (Phi*V).';   
        
end

