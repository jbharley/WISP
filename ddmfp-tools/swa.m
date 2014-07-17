function [V, Phi] = swa( k, d, X, tau, varargin )
%SWA  Sparse wavenumber analysis
%   [DkV Phi] = SWA( K, D, X, TAU, OPTIONS ) recovers the sparse
%   frequency-wavenumber representation of data
%
%   INPUTS:   
%       K: An N-by-1 vector of wavenumbers
%       D: An M-by-1 vector of distancess
%       X: An Q-by-M matrix of frequency-domain data with Q frequencies and 
%          cooresponding to M measurements
%     TAU: For basis pursuit denoising, TAU is the regularization
%          parameter. For orthogonal matching pursuit, TAU is the target
%          sparsity of V. 
%  
%   OPTIONS: 
%           'method': 'bp' for basis pursuit denoising, 'omp' for
%                     orthogonal matching pursuit. Default is 'bp'.
%        'optMatrix': If true, solve the basis pursuit for all frequencies 
%                     at once. Default is false.
%             'plot': If true and optMatrix is false, plots results after 
%                     each frequency. Default is false. 
%
%   OUTPUTS:
%       V: An N-by-Q matrix representing the wave's sparse 
%          frequency-wavenumber representation
%     Phi: An M-by-N matrix representing the medium's propogation matrix
%
%   WARNING: 
%     Use of Basis Pursuit alogrithm requires the CVX toolbox. CVX may
%     downloaded for free here: http://cvxr.com/cvx/.
%   
%   see also: sws, ddmfp, fddmfp
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
    
    % --------------------------------------------------------------------
    % MANAGE INPUT ARGUMENTS
    % --------------------------------------------------------------------

    % CHECK NUMBER OF ARGUMENTS
    if nargin < 4, error('SWA requires 4 or more input arguments.'); end 

    % FIX ARGUMENT FORMATS
    if  iscell(X), X = cell2mat(X); end     % Make a matrix
    
    % SET DEFUALT OPTIONAL ARGUMENTS
    opt.optMatrix     = false;              % Use fully matrix optimization (memory intensive) 
    opt.plot          = false;              % Plot results as they run (doesn't work if optMatrix = true)
    opt.method        = 'bp';               % Optimization method
    
    % PARSE ARGUMENTS AND SIMPLIFY SOME ARGUMENT NAMES
    if ~isempty(varargin), opt = parseArgs(opt, varargin{:}); end
    
    
    % ---------------------------------------------------------------------
    
    % DEFINE LENGTHS
    Q = size(X, 1);  % Number of frequencies
    M = size(d,1);   % Number of measurement
    N = length(k);   % Number of wavenumbers
    
    % NORMALIZE INPUT DATA
    EX = sqrt(sum(abs(X).^2,2));    % Data normalization factors (at each frequency)
    Xn = bsxfun(@times, X, 1./EX);  % Normalized data
    
    % BUILD PROPOGATION MATRICES
    A = exp(-1j*d*k.');                 % Wave propagation frame
    Dr = spdiags(1./sqrt(d), 0, M, M);  % Distance attenuation
    Dk = spdiags(1./sqrt(k), 0, N, N);  % Wavenumber attenutation
    rho = 1./norm(Dr,'fro');            % Normalization: 1./sqrt(sum(1./d))
    Phi = rho*Dr*A;                     % Progagation frame 
                                        %   Note that Phi does not contain 
                                        %   Dk due to column normalization
    
    % COMPUTE SPARSE WAVENUMBER SOLUTION
    switch opt.method
        case 'bp'
            if   opt.optMatrix, DkV = bp(Phi, Xn.', tau);  
            else DkV = cell2mat(arrayfun(@(ii) countLoop(@bp, ii, Q, opt.plot*k, Phi, permute(Xn(ii,:,:), [2 1 3]), tau), 1:Q, 'UniformOutput', false )); end
            
            % DE-BIAS BASIS PURSUIT DENOISING RESULTS
            Y   = (Phi*DkV).';  % Generate denoised signal
            mu  = cell2mat(arrayfun(@(ii) X(ii,:)*Y(ii,:)' ./ norm(Y(ii,:))^2, 1:Q, 'UniformOutput', false ));
            DkV = bsxfun(@times, mu, DkV);  % Frequency-wavenumber representation
            
        case 'omp'
            if   opt.optMatrix, DkV = omp(Phi, X.', tau);  
            else DkV = cell2mat(arrayfun(@(ii) countLoop(@omp, ii, Q, opt.plot*k, Phi, permute(X(ii,:,:), [2 1 3]), tau), 1:Q, 'UniformOutput', false )); end
    end
    
    % REMOVE NORMALIZATION WEIGHTING
    V = Dk\DkV*rho;  
    
    
end


function V = bp( A, X, tau ) 
%BP  Compute the sparse basis pursuit solution
%   V = BP(A, X, TAU) computes the sparse inverse solution V to a basis 
%   matrix A and observations X.
%
%   INPUTS: 
%       A: An M-by-N basis matrix
%       X: An M-by-Q vector of observations
%     tau: Regularization parameter for optimization
%
%  OUTPUTS: 
%       V: The N-by-Q basis pursuit denoising result
%
%   WARNING: 
%     Use of Basis Pursuit alogrithm requires the CVX toolbox. CVX may
%     downloaded for free here: http://cvxr.com/cvx/.
%   
 
    % INITIALIZE VARIABLES
    Q = size(X,2);      % Number of frequencies
    N = size(A,2);      % Number of wavenumbers

    % RUN CVX OPTIMIZATION
    cvx_begin quiet
        variable V(N,Q) complex
        
        % MINIMIZE COST FUNCTION
        minimize( pow_pos(norm( A*V - X, 'fro' ),2) + tau*sum(norms( V, 1, 1 )) );
        
    cvx_end


end



function V = omp( A, X, K ) 
%OMP  Computes the sparse orthogonal matching pursuit solution
%   V = OMP(A, B, K)  computes the sparse solution V given a basis matrix 
%   A and observations X.
%
%   INPUTS: 
%       A: An M-by-N basis matrix
%       X: An M-by-Q vector of observations
%       K: Number of sparse components
%
%  OUTPUTS: 
%       V: The N-by-Q wavenumber basis pursuit denoising results
%

    % DEFINE SIZES 
    N   = size(A,2);            % Number of atoms
    M   = size(A,1);            % Aize of atoms
    Q   = size(X,2);            % Number of outputs/frequencies

    % INITIALIZE VARIABLES
    V         = zeros(N,Q);     % Solution
    V_T       = zeros(K,Q);     % Solution with only non-zero indices
    indx_set  = zeros(K,Q);     % Indices with non-zero values
    atoms     = zeros(M,K,Q);   % Chosen dictionary atoms for each frequency

    % INIATIVE ALGORITHM 
    r   = X;                    % Initial residual
    Vr  = A'*r;                 % Initial solution from residual xr
    
    % LOOP OVER NUMBER OF SPARSE COMPONENTS
    for k = 1:K
        
        % FIND CORRESPONDING INDICES
        [~,ind_new]   = max(abs(Vr), [], 1);             % Find best match
        indx_set(k,:) = ind_new;                         % Add to index set
        atoms(:,k,:)  = permute(A(:,ind_new), [1 3 2]);  % Get cooresponding atom

        % UPDATE RESIDUAL
        for q = 1:Q  % Loop over outputs
            V_T(1:k,q) = atoms(:,1:k,q) \ X(:,q);          % Find least-squares fit
            V( indx_set(1:k,q), q )   = V_T(1:k,q);        % Places results in full vector
            r(:,q)  = X(:,q) - atoms(:,1:k,q)*V_T(1:k,q);  % Find new residual
        end
        
        % COMPUTE SOLUTION FROM RESIDUAL Vr FOR NEXT ITERATION
        if k < K, Vr  = A'*r; end

    end

end

function varargout = countLoop( fnc, n, N, k, varargin )
%COUNTLOOP  Displays the time till complete of a function loop 
%   VARARGOUT = COUNTLOOP( FNC, n, N, k, VARARGIN ) runs the function 
%   FNC, displays the loop count n/N, and displays the estimated time to
%   completion
%
%   INPUTS:   
%       FNC: A function handle to run
%         n: Current interation in loop
%         N: Last interation in loop
%         K: If non-zero, plots the function output with "K" defined as 
%            the horizontal axis
%  VARARGIN: Input parameters to the function "fnc"
%
%  OUTPUTS:
% VARARGOUT: Output of function "fnc"
%


    % DISPLAY TIME INFORMATION
    global tm; if isempty(tm), tm = 0; end
    fprintf('%08i / %08i [Time left: %s]', n, N, datestr(tm/24/3600*(N-n+1), 'HH:MM:SS')); ts = tic; 
    
    % RUN FUNCTION
    M = nargout(fnc); if M == -1, M = 1; end;
    [varargout{1:M}] = fnc(varargin{:});

    % REFRESH TIME INFORMATION
    fprintf(repmat('\b', 1, 41)); 
    
    % REFRESH TIME INFORMATION
    tm = (toc(ts) + tm)/2;
    
    % PLOT RESULT IF REQUESTED
    if any(k)
        if length(varargout{1}) == length(k)
            figure(1);
            plot(k, abs(varargout{1})); axis xy;
            xlabel('Wavenumber [m^{-1}]')
            drawnow;
        elseif length(varargout{1}) == (length(k)^2)
            figure(1)
            V0 = reshape(varargout{1}, length(k), length(k));
            imagesc(k,k,abs(squeeze(V0(:, :))))
            colormap(copper)
            drawnow;
        else
            figure(1);
            plot(abs(varargout{1})); axis xy;
            xlabel('Wavenumber [samples]')
            drawnow;
        end
    end
    
end



function options = parseArgs(options, varargin)

    % ---------------------------------------------------
    % CODE TAKEN FROM STACKOVERFLOW
    %   http://stackoverflow.com/questions/2775263/how-to-deal-with-name-value-pairs-of-function-arguments-in-matlab
    % ---------------------------------------------------
    
    % GET OPTIONS NAMES
    optionNames = fieldnames(options);

    % COUNT ARGUMENTS
    nArgs = length(varargin);
    if round(nArgs/2)~=nArgs/2
       error('SWA needs propertyName/propertyValue pairs')
    end
    
    % PARSE ARGUMENTS
    for pair = reshape(varargin,2,[]) % pair is {propName;propValue}
       inpName = pair{1}; % make case insensitive

       if any(strcmpi(inpName,optionNames))
          % overwrite options. If you want you can test for the right class here
          % Also, if you find out that there is an option you keep getting wrong,
          % you can use "if strcmp(inpName,'problemOption'),testMore,end"-statements
          options.(inpName) = pair{2};
       else
          error('%s is not a recognized parameter name',inpName)
       end
    end

end
