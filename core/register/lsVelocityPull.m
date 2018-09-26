function [ok, v, match, ipsi] = lsVelocityPull(model, dv, r0, match0, mu, f, varargin)
%__________________________________________________________________________
% 
% Performs a line search along a direction to find better initial
% velocity. The line search direction is usually found by Gauss-Newton.
%
%--------------------------------------------------------------------------
%
% FORMAT [ok, ...] = lsVelocityShapePull(model, dv, v0, r0, match0, mu, f, ...)
%
% REQUIRED
% --------
% model  - Structure with fields:
%           * 'name'    : 'normal', 'laplace', 'bernoulli' or 'categorical'
%           * ('sigma2'): Normal variance  [1]
%           * ('b')     : Laplace variance [1]
% dv     - Line search direction (ascent if maximising, descent if minimising)
% v0     - Previous velocity field
% match0 - Previous log-likelihood (matching term)
% mu     - Template (in native space)
% f      - Image (in native space) -> can be a function of `z` (slice number)
%
% KEYWORD ARGUMENTS - FUNCTION HANDLES
% ------------------------------------
% exp      - Function of `v` that recovers the (affine+diffeo) transformation
% pull     - Function of `ipsi` and `mu that warps and reconstruct the template
% llmatch  - Function of `mu` and `f` that computes that matching term
% llreg    - Function of `v` that computes the prior term

% KEYWORD ARGUMENTS - VALUES
% --------------------------
% nit     - Number of line-search iterations [6]
% verbose - Write line search results [false]
% debug   - WWrite debugging stuff [false]
%
% OUTPUT
% ------
% ok    - True if a better parameter value was found
% And if ok == true:
% v     - New initial velocity
% match - New log-likelihood (matching term)
% ipsi  - New complete affine+diffeomorphic mapping
%
%--------------------------------------------------------------------------
% 
% The objective function is supposed to be of the form
%   E = logMatch(f|mu,v) + logN(v | mean, inv(lam*L))
%
% The regularisation is split in two parts (lam and L), because L (alone)
% is used to exponentiate the velocity field (i.e., it also
% parameterises shooting.)
%__________________________________________________________________________

    % ---------------------------------------------------------------------
    % Parse inputs
    % ------------
    
    p = inputParser;
    p.FunctionName = 'lsVelocity';
    p.addRequired('model',  @(X) ischar(X) || (isstruct(X) && isfield(X, 'name')));
    p.addRequired('dv',     @checkarray);
    p.addRequired('v0',     @checkarray);
    p.addRequired('match0', @isscalar);
    p.addRequired('mu',     @checkarray);
    p.addRequired('f',      @checkarray);
    p.addParameter('exp',     @(X) isa(X, 'function_handle'), @(v) exponentiateVelocity(v, 'iphi'));
    p.addParameter('pull',    @(X) isa(X, 'function_handle'), @pullTemplate);
    p.addParameter('llmatch', @(X) isa(X, 'function_handle'), @(mu,f) llMatching('normal', mu, f));
    p.addParameter('llregr',  @(X) isa(X, 'function_handle'), @(r) llPriorVelocity(r, 'fast'));
    p.addParameter('llregv',  @(X) isa(X, 'function_handle'), @(v) llPriorVelocity(v, 'fast'));
    p.addParameter('v0',   [],      @checkarray);
    p.addParameter('reg0', nan,     @isscalar);
    p.addParameter('nit',  6,     @isscalar);
    p.addParameter('verbose', false, @isscalar);
    p.addParameter('debug',   false, @isscalar);
    
    p.parse(model, dv, r0, match0, mu, f, varargin{:});
    
    exponentiate = p.Results.exp;
    pull         = p.Results.pull;
    llmatch      = p.Results.llmatch;
    llreg        = p.Results.llreg;
    v0           = p.Results.v0;
    reg0         = p.Results.reg0;
    nit          = p.Results.nit;
    verbose      = p.Results.verbose;
    debug        = p.Results.debug;
    % ---------------------------------------------------------------------
    
    if debug, fprintf('* lsVelocityShapePull\n'); end
    
    % ---------------------------------------------------------------------
    % Load some data (in case it is on disk)
    v0 = numeric(v0);
    dv = numeric(dv);
    
    % ---------------------------------------------------------------------
    % Set some default parameter value
    if isnan(vel0), reg0 = llreg(v0);
    end
    
    % ---------------------------------------------------------------------
    armijo = 1;             % Armijo factor
    ll0    = match0 + reg0; % Log-likelihood (only parts that depends on v)
    if verbose
        printInfo('header');
        printInfo('initial', ll0, match0, reg0);
    end
    
    % ---------------------------------------------------------------------
    % Loop
    for i=1:nit
        % -----------------------------------------------------------------
        % Compute transform
        r    = r0 + dv / armijo;
        v    = v0 + dv / armijo;
        ipsi = exponentiate(v);
        if nargout < 4
            clear ipsi
        end
        
        % -----------------------------------------------------------------
        % Matching term
        match = 0;
        for z=1:size(ipsi,3)
            wmu = pull(ipsi(:,:,z,:), mu);
            if isa(f, 'function_handle'), fz = f(z);
            else,                         fz = f(:,:,z,:);
            end
            match = match + llmatch(wmu, fz);
            clear wmu fz
        end
        
        % -----------------------------------------------------------------
        % Prior term
        reg = llreg(v);
        
        % -----------------------------------------------------------------
        % Complete objective function
        ll = match + reg;
        
        if verbose, printInfo(i, ll0, match, reg); end
        
        if ll <= ll0
            if verbose, printInfo('failed'); end
            armijo = armijo * 2;
        else
            if verbose, printInfo('success'); end
            ok = true;
            return
        end
    end
    
    % ---------------------------------------------------------------------
    % Complete failure
    if verbose, printInfo('end'); end
    ok = false;

end

function printInfo(which, oll, llm, llr)
    if ischar(which) 
        if strcmpi(which, 'header')
            fprintf('\nV - LineSearch | Armijo  | %12s = %12s + %12s | %12s\n', 'RLL', 'Match', 'Reg', 'Diff');
        elseif strcmpi(which, 'initial')
            fprintf('V - LineSearch | Initial | %12.4f = %12.4f + %12.4f \n', oll, llm, llr);
        elseif strcmpi(which, 'failed')
            fprintf('| Failed\n');
        elseif strcmpi(which, 'success')
            fprintf('| Success\n\n');
            fprintf(repmat(' ',1,50));
        elseif strcmpi(which, 'end')
            fprintf('V - LineSearch | Complete failure\n\n');
            fprintf(rempat(' ',1,50));
        end
    else
        fprintf('V - LineSearch | Try %3d | %12.4f = %12.4f + %12.4ff | %12.4f ', which, llm+llr, llm, llr, llm+llr-oll);
    end
end
        