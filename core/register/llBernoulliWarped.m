function ll = llBernoulliWarped(mu, f, varargin)
% FORMAT ll = llBernoulliWarped(mu, f, (type), ('loop', loop), ('par', par))
%
% ** Required **
% mu   - Template warped in image space.
% f    - Observed image in native space.
% ** Optional **
% type - Type of the provided template:
%        * 'proba': Probability space [default]
%        * 'log'  : Log-probability space
% ** Keyword arguments **
% loop - Specify how to split data processing
%        ('slice' or 'none' [default])
% par  - If true, parallelise processing [default: false]
% ** Output **
% ll   - Log-likelihood of the Normal matching term

    % --- Parse inputs
    p = inputParser;
    p.FunctionName = 'llBernoulliWarped';
    p.addRequired('mu',   @checkarray);
    p.addRequired('f',    @checkarray);
    p.addOptional('type',    'proba', @(X) ischar(X) && any(strcmpi(X, {'log', 'proba'})) );
    p.addParameter('loop',   '',      @ischar);
    p.addParameter('par',    false,   @isscalar);
    p.addParameter('debug',  false,   @isscalar);
    p.parse(mu, f, varargin{:});
    type  = p.Results.type;
    par   = p.Results.par;
    loop  = p.Results.loop;
    
    if p.Results.debug, fprintf('* llBernoulliWarped\n'); end;
    
    % --- Optimise parallelisation and splitting schemes
    [par, loop] = autoParLoop(par, loop, isa(mu, 'file_array'), size(mu, 3), 1);
    
    % --- Read dimensions
    dim  = [size(mu) 1 1];
    dlat = dim(1:3);
    
    % --- Initialise
    ll = 0;
    
    % --- No loop
    if strcmpi(loop, 'none')
        if p.Results.debug, fprintf('   - No loop\n'); end;
        if strcmpi(type, 'proba')
            ll = ll + onMemoryProba(mu(:,:,:,1), f(:,:,:,1));
        elseif strcmpi(type, 'log')
            ll = ll + onMemoryLog(mu(:,:,:,1), f(:,:,:,1));
        else
            error('Unknown template type %s', type);
        end
        
    % --- Loop on slices
    elseif strcmpi(loop, 'slice')
        if p.Results.debug
            if par > 0
                fprintf('   - Parallelise on slices\n'); 
            else
                fprintf('   - Serialise on slices\n'); 
            end
        end
        if strcmpi(type, 'proba')
            parfor (z=1:dlat(3), par)
                ll = ll + onMemoryProba(mu(:,:,z,1), f(:,:,z,1));
            end
        elseif strcmpi(type, 'log')
            parfor (z=1:dlat(3), par)
                ll = ll + onMemoryLog(mu(:,:,z,1), f(:,:,z,1));
            end
        else
            error('Unknown template type %s', type);
        end
        
    end
end

function ll = onMemoryProba(mu, f)
    mu    = single(numeric(mu));
    f     = single(numeric(f));
    msk   = isfinite(f) & isfinite(mu);
    mu    = max(eps('single'), min(1 - eps('single'), mu));
    ll    = sumall(f(msk) .* log(mu(msk)) + (1 - f(msk)) .* log(1 - mu(msk)));
end

function ll = onMemoryLog(a, f)
    a     = single(numeric(a));
    f     = single(numeric(f));
    msk   = isfinite(f) & isfinite(a);
    ll    = sumall(f(msk) .* a(msk) - log(1 + exp(a(msk))));
end