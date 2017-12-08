function ll = llNormalWarped(mu, f, varargin)
% FORMAT ll = llNormalWarped(mu, f, (s), ('loop', loop), ('par', par))
%
% ** Required **
% mu   - Template warped in image space.
% f    - Observed image in native space.
% ** Optional **
% s    - Noise variance (one per modality).
% ** Keyword arguments **
% loop - Specify how to split data processing
%        ('component', 'slice' or 'none' [default])
% par  - If true, parallelise processing [default: false]
% ** Output **
% ll   - Log-likelihood of the Normal matching term

    % --- Parse inputs
    p = inputParser;
    p.FunctionName = 'llNormalWarped';
    p.addRequired('mu',  @checkarray);
    p.addRequired('f',   @checkarray);
    p.addOptional('s',   1,  @(X) isscalar(X) || checkarray(X));
    p.addParameter('loop',   '',    @ischar);
    p.addParameter('par',    false, @isscalar);
    p.addParameter('debug',  false, @isscalar);
    p.parse(mu, f, varargin{:});
    par = p.Results.par;
    loop = p.Results.loop;
    
    if p.Results.debug, fprintf('* llNormalWarped\n'); end;
    
    % --- Optimise parallelisation and splitting schemes
    [par, loop] = autoParLoop(par, loop, isa(mu, 'file_array'), size(mu, 3), size(mu, 4));
    
    % --- Prepare sigma
    nc = size(mu, 4);
    s = p.Results.s;
    if length(s) == 1
        s = s * ones(1, nc);
    end
    
    % --- Read dimensions
    dim  = [size(mu) 1 1];
    dlat = dim(1:3);
    
    % --- Initialise
    ll = 0;
    
    % --- No loop
    if strcmpi(loop, 'none')
        if p.Results.debug, fprintf('   - No loop\n'); end;
        ll = ll + onMemory(mu, f, s);
    
    % --- Loop on components
    elseif strcmpi(loop, 'component')
        if p.Results.debug
            if par > 0
                fprintf('   - Parallelise on components\n'); 
            else
                fprintf('   - Serialise on components\n'); 
            end
        end
        if ~par
            for k=1:nc
                ll = ll + onMemory(mu(:,:,:,k), f(:,:,:,k), s(k));
            end
        elseif isa(mu, 'file_array') && isa(f, 'file_array')
            parfor (k=1:nc, par)
                ll = ll + onMemory(slicevol(mu, k, 4), slicevol(f, k, 4), s(k));
            end
        elseif isa(mu, 'file_array')
            parfor (k=1:nc, par)
                ll = ll + onMemory(slicevol(mu, k, 4), f(:,:,:,k), s(k));
            end
        elseif isa(f, 'file_array')
            parfor (k=1:nc, par)
                ll = ll + onMemory(mu(:,:,:,k), slicevol(f, k, 4), s(k));
            end
        else
            parfor (k=1:nc, par)
                ll = ll + onMemory(mu(:,:,:,k), f(:,:,:,k), s(k));
            end
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
        if ~par
            for z=1:dlat(3)
                ll = ll + onMemory(mu(:,:,z,:), f(:,:,z,:), s);
            end
        elseif isa(mu, 'file_array') && isa(f, 'file_array')
            parfor (z=1:dlat(3), par)
                ll = ll + onMemory(slicevol(mu, z, 3), slicevol(f, z, 3), s);
            end
        elseif isa(mu, 'file_array')
            parfor (z=1:dlat(3), par)
                ll = ll + onMemory(slicevol(mu, z, 3), f(:,:,z,:), s);
            end
        elseif isa(f, 'file_array')
            parfor (z=1:dlat(3), par)
                ll = ll + onMemory(mu(:,:,z,:), slicevol(f, z, 3), s);
            end
        else
            parfor (z=1:dlat(3), par)
                ll = ll + onMemory(mu(:,:,z,:), f(:,:,z,:), s);
            end
        end
        
    end
end

function ll = onMemory(mu, f, s)
    nc    = size(mu, 4);
    mu    = single(numeric(mu));
    f     = single(numeric(f));
    msk   = isfinite(f) & isfinite(mu);
    count = sum(sum(sum(msk, 3), 2), 1);
    ll    = 0;
    for k=1:nc
        ll = ll - 0.5 * count(k) * (s(k) + log(2 * pi)) ...
                - 0.5 / s(k) * sumall((mu(msk(:,:,:,k)) - f(msk(:,:,:,k))).^2);
    end
end
