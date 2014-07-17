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
% Code: Written by: Joel B. Harley
% Last updated: July 16, 2014
% -------------------------------------------------------------------------
% If this code is used for a research publication, please cite:
% J.B. Harley, J.M.F. Moura, "Sparse Recovery of the Multimodal and 
% Dispersive Characteristics of Lamb Waves," Journal of the Acoustical
% Society of America, vol. 133, no. 5, May 2013.
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
    A = exp(-1j*d*k.');                 % Complex exponentials
    Dr = spdiags(1./sqrt(d), 0, M, M);  % Distance attenuation
    Dk = spdiags(1./sqrt(k), 0, N, N);  % Wavenumber attenutation 
    Phi = Dr*A*Dk;                      % Progagation frame 
                                        
    % SYNTHESIZE SIGNALS
    X = (Phi*V).';   
        
end

