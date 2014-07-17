function [pxl, label] = ddmfp( k, d, fn, x, V )
% DDMFP  Data-driven matched field processing
%   PXL = DDMFP( K, D FN, X, V, OPTIONS ) performs data-driven matched  
%   field processing to localize a target
%
%   INPUTS: 
%       K: An N-by-1 vector of wavenumbers, corresponding to the  
%          wavenumbers represented by parameter V
%       D: An M-by-L matrix of distances associated with L grid points and
%          M measurements
%      FN: A Qn-by-1 matrix of frequencies to perform localization over
%       X: A Q-by-M matrix of time-domain signals with Q samples and 
%          cooresponding to M measurements
%       V: An N-by-Q matrix representing a signal's sparse 
%          frequency-wavenumber representation
%
%   OPTIONS: 
%       'incoherent': If true, DDMFP uses the incoherent matched field
%                     processor (true by default)
%         'coherent': If true, DDMFP uses the coherent matched field
%                     processor (true by default)
%
%   OUTPUTS:
%     PXL: An L-by-P matrix of P ambiguity sufaces, corresponding to P 
%          processors
%
%   see also: swa, sws, fddmfp, mfp, fmfp
%

% -------------------------------------------------------------------------
% Code written by: Joel B. Harley
% Last updated: July 16, 2014
% -------------------------------------------------------------------------
% If this code is used for a research publication, please cite:
% J.B. Harley, J.M.F. Moura, "Data-driven matched field processing for Lamb
% wave structural health monitoring," Journal of the Acoustical Society of
% America, vol. 135, no. 3, March 2014.
% -------------------------------------------------------------------------
%

    % CHECK NUMBER OF ARGUMENTS
    if nargin < 4, error('DDMFP requires 4 or more input arguments.'); end 
    
    % FORCE VECTORS TO BE COLUMN VECTORS
    fn = fn(:);
    k  = k(:);
    
    % SET DEFUALT OPTIONAL ARGUMENTS
    opt.incoherent = true;
    opt.coherent   = true;
    
    % PARSE ARGUMENTS 
    if ~isempty(varargin), opt = parseArgs(opt, varargin{:}); end 
    
    % CHECK FOR ERRORS IN PARAMETERS
    if size(d,1) ~= size(x,2), error('Check dimensions of D and X'); end
    
    
    % ---------------------------------------------------------------------
    
    % INITIALIZE VARIABLES
    L = size(d,2);       % Number of grid points
    M = size(d,1);       % Number of distances (sensor pairs)
    Qn = size(fn,1);     % Number of time samples in data
    P = opt.incoherent + opt.coherent;  % Number of processors
    
    % INITIALIZE RESULTS
    pxl   = zeros(L,P);  % Pixel vector
    label = cell(P,1);   % Processor labels
    
    % COMPUTE FOURIER TRANSFORM OF DATA
    X = fft(x).';        % Fourier transform
    X = X(:,fn);         % Use selected frequencies
    
    % CONSTRUCT CROSS-SPECTRAL DENSITY MATRIX AND ITS INVERSE
    KK = zeros(M,M,Qn); KKi = zeros(M,M,Qn);
    for n = 1:Qn
        KK(:,:,n) = X(:,n)*X(:,n)' + opt.loadingFactor*eye(M);
        KKi(:,:,n) = pinv(KK(:,:,n));
    end
    
    % PERFORM DATA-DRIVEN MATCHED FIELD PROCESSING AT EACH GRID POINT
    tm = 0; fprintf(repmat(' ', 1, 41));
    for n = 1:L 
        fprintf([repmat('\b', 1, 41) '%08i / %08i [Time left: %s]'], n, L, datestr(tm/24/3600*(L-n+1), 'HH:MM:SS')); ts = tic; 
        
        % PREDICT AND NORMALIZE DATA
        Y = sws( k, d(:,n), V ).';   % Predicted signal
        nm = 1./sqrt(sum(abs(Y).^2)); nm(isinf(nm)) = 0;
        Yn = bsxfun(@times, Y, nm);  % Normalized signal
        
        % DEFINE LOCALIZATION PROCESSORS
        p = 1;
        if opt.incoherent, pxl(n,p) = real(sum(diag(abs(Yn'*X).^2))); label{p} = 'Incoherent'; p=p+1; end  % Incoherent Matched Field Processor
        if opt.coherent,   pxl(n,p) = abs(trace(Y'*X)).^2./(norm(Y, 'fro')).^2; label{p} = 'Coherent'; p=p+1; end  % Coherent Matched Field Processor
        
        % REFRESH TIME INFORMATION
        tm = (toc(ts) + tm)/2;
    end
    fprintf(repmat('\b', 1, 41));

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
       error('DDMFP needs propertyName/propertyValue pairs')
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
