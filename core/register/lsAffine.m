function [ok, q, llm, llq, A, pf, c, ipsi] = lsAffine(model, dq, q0, llm0, mu, f, varargin)
% FORMAT [ok, q, llm, llq, A, pf, c, ipsi] = lsAffine(model, dq, q0, llm0, mu, f, ...)
%
% ** Required **
% model - Structure with fields:
%           * 'name'    : 'normal', 'laplace', 'bernoulli' or 'categorical'
%           * ('sigma2'): Normal variance  [1]
%           * ('b')     : Laplace variance [1]
% dq   - Line search direction (ascent if maximising, descent if minimising)
% q0   - Previous parameter value
% llm0 - Previous log-likelihood (matching term)
% mu   - Template (in native space)
% f    - Image (in native space)
% ** Keyword arguments **
% llq0 - Previous log-likelihood (prior term) [compute]
% B    - Affine basis [affine_basis('affine')]
% regq - Precision matrix of the affine parameters [none]
% rind - Indices of regularised affine parameters []
% iphi - Diffeomorphic part [identity]
% Mf   - Image voxel-to-world mappinf [eye(4)]
% Mmu  - Template voxel-ti-world mapping [eye(4)]
% nit  - Number of line-search iterations [6]
% loop - How to split processing [auto]
% par  - If true, parallelise processing [false]
% ** Output **
% ok   - True if a better parameter value was found
% q    - New parameter value
% llm  - New log-likelihood (matching term)
% llq  - New log-likelihood (prior term)
% A    - New affine transform
% pf   - New pushed image
% c    - New pushed voxel count
% ipsi - (if needed) New complete affine+diffeomorphic mapping
%
% Performs a line search along a direction to find better affine
% parameters. The line search direction is usually found by Gauss-Newton

    % --- Parse inputs
    p = inputParser;
    p.FunctionName = 'lsAffine';
    p.addRequired('model',  @(X) isstruct(X) && isfield(X, 'name'));
    p.addRequired('dq',     @checkarray);
    p.addRequired('q0',     @checkarray);
    p.addRequired('llm0',   @isscalar);
    p.addRequired('mu',     @checkarray);
    p.addRequired('f',      @checkarray);
    p.addParameter('llq0',     nan,    @iscalar);
    p.addParameter('B',        [],     @checkarray);
    p.addParameter('rind',     [],     @isnumeric);
    p.addParameter('regq',     []);
    p.addParameter('iphi',     [],     @checkarray);
    p.addParameter('Mf',       eye(4), @(X) isnumeric(X) && issame(size(X), [4 4]));
    p.addParameter('Mmu',      eye(4), @(X) isnumeric(X) && issame(size(X), [4 4]));
    p.addParameter('nit',      6,      @isscalar);
    p.addParameter('par',      false,  @isscalar);
    p.addParameter('loop',     '',     @(X) ischar(X) && any(strcmpi(X, {'slice', 'component', 'none', ''})));
    p.addParameter('output',   []);
    p.addParameter('verbose',  false,  @isscalar);
    p.addParameter('debug',    false,  @isscalar);
    p.parse(model, dq, q0, llm0, mu, f, varargin{:});
    llq0   = p.Results.llq0;
    B      = p.Results.B;
    rind   = p.Results.rind;
    regq   = p.Results.regq;
    iphi   = p.Results.iphi;
    Mf     = p.Results.Mf;
    Mmu    = p.Results.Mmu;
    nit    = p.Results.nit;
    par    = p.Results.par;
    loop   = p.Results.loop;
    output = p.Results.output;
    verbose = p.Results.verbose;
    debug  = p.Results.debug;
    
    if debug, fprintf('* lsAffine\n'); end;
    
    % --- Set some default parameter value
    if isempty(B)
        error('Missing affine basis')
    end
    if ~isempty(regq) && isempty(rind)
        rind = (1:size(regq,1)) + numel(q0) - size(regq, 1);
    end
    if isnan(llq0)
        if ~isempty(regq),  llq0 = llPriorAffine(q0(rind), regq, 'fast', 'debug', debug);
        else                llq0 = 0; end;
    end
    dim  = [size(mu) 1 1];
    lat = dim(1:3);
    if isempty(iphi)
        iphi = warps('identity', lat);
    end
    
    % --- Load some data (in case it is on disk)
    q0 = numeric(q0);
    dq = numeric(dq);
    
    % --- Initialise line search
    armijo = 1;           % Armijo factor
    ok     = false;       % Found a better ll ?
    ll0    = llm0 + llq0; % Log-likelihood (only parts that depends on q)
    dimf   = [size(f) 1 1];
    latf   = dimf(1:3);
    
    if verbose
        printInfo('header');
        printInfo('initial', ll0, llm0, llq0);
    end
    
    % //!\\ Here, the pushed image is temporarily kept in memory. This
    % might use up a large ammount of RAM when a lot of classes are used.
    % I might need to store it in a temporary file_array.
    
    % --- Loop
    for i=1:nit
        q = q0 + dq / armijo;
        A = exponentiateAffine(q, B, 'debug', debug);
        ipsi = reconstructIPsi(A, iphi, 'lat', latf, 'Mf', Mf, 'Mmu', Mmu, 'debug', debug);
        [pf, c] = pushImage(ipsi, f, lat, 'par', par, 'loop', loop, 'debug', debug);
        llm = llMatching(model, mu, pf, c, 'par', par, 'loop', loop, 'debug', debug);
        if checkarray(regq)
            llq = llPriorAffine(q(rind), regq, 'fast', 'debug', debug);
        else
            llq = 0;
        end
        ll  = llm + llq;
        
        if verbose, printInfo(i, ll0, llm, llq); end;
        
        if ll <= ll0
            if verbose, printInfo('failed'); end;
            armijo = armijo * 2;
        else
            if verbose, printInfo('success'); end;
            if checkarray(regq)
                llq = llPriorAffine(q(rind), regq, 'debug', debug);
            end
            ok  = true;
            return
        end
    end
    
    if verbose, printInfo('end'); end;
    q   = q0;

end

function printInfo(which, oll, llm, llq)
    if ischar(which) 
        if strcmpi(which, 'header')
            fprintf('Q - LineSearch | Armijo  | %12s = %12s + %12s | %12s\n', 'RLL', 'LL-Match', 'RLL-Prior', 'LL-Diff');
        elseif strcmpi(which, 'initial')
            fprintf('Q - LineSearch | Initial | %12.6f = %12.6f + %12.6f \n', oll, llm, llq);
        elseif strcmpi(which, 'failed')
            fprintf('| Failed\n');
        elseif strcmpi(which, 'success')
            fprintf('| Success\n');
        elseif strcmpi(which, 'end')
            fprintf('Q - LineSearch | Complete failure\n');
        end
    else
        fprintf('Q - LineSearch | Try %3d | %12.6f = %12.6f + %12.6f | %12.6f ', which, llm+llq, llm, llq, llm+llq-oll);
    end
end
        