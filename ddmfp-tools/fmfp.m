function [pxl, label] = fmfp( c, d, fn, x, s, varargin )
% FMFP  Fast matched field processing (delay-and-sum localization)
%   PXL = FMFP( C, D, FN, X ) perform delay-sum matched field localization
%   with velocity C, grid distances D, frequencies FN, data X, and
%   exicitation S
%
%   PXL = FMFP( C, D, FN, X, 'envelope-delay' ) perform delay-sum matched 
%   field localization with enveloped signals
%
%   INPUTS: 
%       C: Group velocity of signal of interest (should be normalized by 
%          the sampling rate so that is has units [m / sample])
%       D: An M-by-L matrix of distances associated with L grid points and
%          M measurements
%      FN: A Qn-by-1 matrix of frequencies to perform localization over
%       X: A Q-by-M matrix of time-domain signals with Q samples and 
%          cooresponding to M measurements
%       S: A Q-by-1 vector of the time-domain excitation signal
%
%   OUTPUTS:
%     PXL: An L-by-1 vector corresponding to an ambiguity suface
%
%   WARNING: 
%     This "fast" version of MFP may be memory intensive for a large number
%     of distances M in D and/or large number of samples Q in X
%   
%   see also: mfp, ddmfp, fddmfp
%

% -------------------------------------------------------------------------
% Code written by: Joel B. Harley
% Last updated: July 16, 2014
% -------------------------------------------------------------------------
%

    % --------------------------------------------------------------------
    % MANAGE INPUT ARGUMENTS
    % --------------------------------------------------------------------

    % CHECK NUMBER OF ARGUMENTS
    if nargin < 4, error('SWA requires 4 or more input arguments.'); end 

    % FIX ARGUMENT FORMATS
    if  iscell(x), x = cell2mat(x); end     % Make a matrix
    c = c(:);                               % Make a column vector
    fn = fn(:);                             % Make a column vector
    
    % SET DEFUALT OPTIONAL ARGUMENTS
    opt.loadingFactor = 0.1;
    opt.model = 'delay';
    opt.coherent = true;
    opt.beta = 1;
    
    % PARSE ARGUMENTS 
    if ~isempty(varargin), opt = parseArgs(opt, varargin{:}); end 


    % ---------------------------------------------------------------------   
    
    % INITIALIZE VARIABLES
    L = size(d,2);       % Number of grid points
    M = size(d,1);       % Number of distances (sensor pairs)
    Q = size(x,1);       % Number of time samples in data
    P = opt.coherent;    % Number of processors
    p = 1;
    
    % INITIALIZE RESULTS
    pxl   = zeros(L,P);  % Pixel vector
    label = cell(P,1);   % Processor labels
    
    % BUILD RANGE
    Lm = L*opt.beta;                % Number of ranges 
                                    % (times beta to is arbitrarily used to improve accuracy)
    Dm = max(max(d));               % Maximum distance
    r = linspace(Dm/Lm, Dm, Lm).';  % Ranges 
    
    % COMPUTE FREQUENCY AXIS
    fc = floor(Q/2)+1; f = ifftshift((((1:Q)-fc)/(Q)));
    
    % COMPUTE FOURIER TRANSFORM
    S = fft(s, Q);
    if strcmp(opt.model, 'envelope-delay')
        X = fft(abs(hilbert(x))).';
    else
        X = fft(x).';
    end
    X = X(:,fn); 
    
    % DEFINE MODEL
    if strcmp(opt.model, 'delay')
        Y = bsxfun(@times, S(fn).', bsxfun(@times, 1./sqrt(r), exp(-1j*2*pi*r/c*f(fn))));  % Build wave frame
    elseif strcmp(opt.model, 'envelope-delay')
        Y0 = bsxfun(@times, S(1:Q).', bsxfun(@times, 1./sqrt(r), exp(-1j*2*pi*r/c*f(1:Q))));  % Build wave frame
        Y = fft(abs(hilbert(real(ifft(Y0.'))))).'; Y = Y(:,fn); 
    end
    
    % --------------------------------------------------------------------
    % COHERENT MATCHED FIELD PROCESSOR
    % --------------------------------------------------------------------
    if opt.coherent
        % LOOP OVER NUMBER OF SENSOR PAIRS
        N0 = zeros(L,1); D0 = zeros(L,1);
        tm = 0; fprintf(repmat(' ', 1, 41));
        for m = 1:M
           fprintf([ repmat('\b', 1, 41) '%08i / %08i [Time left: %s]'], m, M, datestr(tm/24/3600*(M-m+1), 'HH:MM:SS')); ts = tic;    
           N0 = N0 + interp1(r,Y*X(m,:)'       ,d(m,:)).';  % numerator
           D0 = D0 + interp1(r,sum(abs(Y).^2,2),d(m,:)).';  % denominator
           tm = (toc(ts) + tm)/2;
        end
        fprintf(repmat('\b', 1, 41));
        pxl(:,p) = abs(N0).^2 ./ D0;  % Coherent matched field processor
        label{p} = 'Coherent'; p=p+1;
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

