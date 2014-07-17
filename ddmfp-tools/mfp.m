function pxl = mfp( c, d, fn, x, s, varargin )
% MFP  Matched field processing (delay-and-sum localization)
%   PXL = MFP( C, D, FN, X ) perform delay-sum matched field localization
%   with velocity C, grid distance D, frequencies FN, data X, and
%   exicitation S
%
%   PXL = MFP( C, D, FN, X, 'envelope-delay' ) perform delay-sum matched 
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
%   see also: fmfp, ddmfp, fddmfp
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
    opt.model = 'delay';
    
    % PARSE ARGUMENTS 
    if ~isempty(varargin), opt = parseArgs(opt, varargin{:}); end 


    % --------------------------------------------------------------------
    % COMPUTE SIGNAL AND CROSS SPECTRAL DENSITY MATRIx
    % --------------------------------------------------------------------    
    
    % COMPUTE FOURIER TRANSFORM
    if strcmp(opt.model, 'envelope-delay')
        X = fft(abs(hilbert(x))).';
    else
        X = fft(x).';
    end
    Q = size(X,2);
    X = X(:,fn);
    S = fft(s, Q);
    
    % FREQUENCY AXIS
    r = floor(Q/2)+1; f = ifftshift((((1:Q)-r)/(Q)));
    
    % --------------------------------------------------------------------
    % PERFORM FUNCTION
    % --------------------------------------------------------------------    
    
    % INITIALIZE VARIABLES
    L  = size(d,2);    % Number of grid points
    pxl = zeros(L,1);  % Pixel vector
    tm  = 0;           % Initialize time
    
    % PERFORM LOCALIZATION AT EACH GRID POINT
    fprintf(repmat(' ', 1, 41));
    for n = 1:L 
        fprintf([ repmat('\b', 1, 41) '%08i / %08i [Time left: %s]'], n, L, datestr(tm/24/3600*(L-n+1), 'HH:MM:SS')); ts = tic; 
        
        % DEFINE MODEL
        if strcmp(opt.model, 'delay')
            Y = bsxfun(@times, S(fn).', bsxfun(@times, 1./sqrt(d(:,n)), exp(-1j*2*pi*d(:,n)/c*f(fn))));  % Build wave frame
        elseif strcmp(opt.model, 'envelope-delay')
            Y0 = bsxfun(@times, S(1:Q).', bsxfun(@times, 1./sqrt(d(:,n)), exp(-1j*2*pi*d(:,n)/c*f(1:Q))));  % Build wave frame
            Ye = fft(abs(hilbert(real(ifft(Y0.'))))).'; Y = Ye(:,fn); 
        end

        % DEFINE LOCALIZATION PROCESSORS
        pxl(n,1) = real(trace(Y'*X))./(norm(Y, 'fro'));                                % "Coherent" Matched Field Processor

        % REFRESH TIME INFORMATION
        tm = (toc(ts) + tm)/2;
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
