function ll = llBernoulliPushed(mu, f, c, varargin)
%__________________________________________________________________________
%
% Log-likelihood of the Bernoulli matching term
%
% -------------------------------------------------------------------------
%
% FORMAT ll = llBernoulliPushed(mu, f, c, (type), ...)
%
% REQUIRED
% --------
% mu   - Template
% f    - Observed image pushed in the template space.
% c    - Pushed voxel count.
% 
% OPTIONAL
% --------
% type - Type of the provided template:
%        * 'proba': Probability space [default]
%        * 'log'  : Log-probability space
%
% KEYWORD ARGUMENTS
% -----------------
% bb   - Bounding box (if different between template and pushed image)
% loop - Specify how to split data processing
%        ('component', 'slice' or 'none' [default])
% par  - If true, parallelise processing [default: false]
% 
% OUTPUT
% ------
% ll   - Log-likelihood of the Bernoulli matching term
%__________________________________________________________________________

    % --- Parse inputs
    p = inputParser;
    p.FunctionName = 'llBernoulliPushed';
    p.addRequired('mu',   @checkarray);
    p.addRequired('f',    @checkarray);
    p.addRequired('c',    @checkarray);
    p.addOptional('type',    'proba', @(X) ischar(X) && any(strcmpi(X, {'log', 'proba'})) );
    p.addParameter('bb',     struct, @isstruct);
    p.addParameter('loop',   '',      @ischar);
    p.addParameter('par',    false,   @isscalar);
    p.addParameter('debug',  false,   @isscalar);
    p.parse(mu, f, c, varargin{:});
    type  = p.Results.type;
    bb    = p.Results.bb;
    par   = p.Results.par;
    loop  = p.Results.loop;
    
    if p.Results.debug, fprintf('* llBernoulliPushed\n'); end;
    
    % --- Optimise parallelisation and splitting schemes
    [par, loop] = autoParLoop(par, loop, isa(mu, 'file_array'), size(mu, 3), size(mu, 4));
    
    % --- Default bounding box
    if ~isfield(bb, 'x')
        bb.x = 1:size(mu, 1);
    end
    if ~isfield(bb, 'y')
        bb.y = 1:size(mu, 2);
    end
    if ~isfield(bb, 'z')
        bb.z = 1:size(mu, 3);
    end
    
    % --- Read dimensions
    dim  = [numel(bb.x)  numel(bb.y) numel(bb.z)];
    dlat = dim(1:3);
    
    % --- Initialise
    ll = 0;
    
    % --- No loop
    if strcmpi(loop, 'none')
        if p.Results.debug, fprintf('   - No loop\n'); end;
        if strcmpi(type, 'proba')
            ll = ll + onMemoryProba(mu, f, c);
        elseif strcmpi(type, 'log')
            ll = ll + onMemoryLog(mu, f, c);
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
            if strcmpi(type, 'proba')
                for z=1:dlat(3)
                    ll = ll + onMemoryProba(mu(bb.x,bb.y,z,:), f(:,:,z,:), c(:,:,z,:));
                end
            elseif strcmpi(type, 'log')
                for z=1:dlat(3)
                    ll = ll + onMemoryLog(mu(bb.x,bb.y,z,:), f(:,:,z,:), c(:,:,z,:));
                end
            end
        elseif isa(mu, 'file_array') && isa(f, 'file_array')
            if strcmpi(type, 'proba')
                parfor (z=1:dlat(3), par)
                    ll = ll + onMemoryProba(slicevol(mu, {bb.x, bb.y, z}), slicevol(f, z, 3), slicevol(c, z, 3), b);
                end
            elseif strcmpi(type, 'log')
                parfor (z=1:dlat(3), par)
                    ll = ll + onMemoryLog(slicevol(mu, {bb.x, bb.y, z}), slicevol(f, z, 3), slicevol(c, z, 3), b);
                end
            end
        elseif isa(mu, 'file_array')
            if strcmpi(type, 'proba')
                parfor (z=1:dlat(3), par)
                    ll = ll + onMemoryProba(slicevol(mu, {bb.x, bb.y, z}), f(:,:,z,:), c(:,:,z,:));
                end
            elseif strcmpi(type, 'log')
                parfor (z=1:dlat(3), par)
                    ll = ll + onMemoryLog(slicevol(mu, {bb.x, bb.y, z}), f(:,:,z,:), c(:,:,z,:));
                end
            end
        elseif isa(f, 'file_array')
            mu = mu(bb.x, bb.y, bb.z,:);
            if strcmpi(type, 'proba')
                parfor (z=1:dlat(3), par)
                    ll = ll + onMemoryProba(mu(:,:,z,:), slicevol(f, z, 3), slicevol(c, z, 3));
                end
            elseif strcmpi(type, 'log')
                parfor (z=1:dlat(3), par)
                    ll = ll + onMemoryLog(mu(:,:,z,:), slicevol(f, z, 3), slicevol(c, z, 3));
                end
            end
        else
            mu = mu(bb.x, bb.y, bb.z,:);
            if strcmpi(type, 'proba')
                parfor (z=1:dlat(3), par)
                    ll = ll + onMemoryProba(mu(:,:,z,:), f(:,:,z,:), c(:,:,z,:));
                end
            elseif strcmpi(type, 'log')
                parfor (z=1:dlat(3), par)
                    ll = ll + onMemoryLog(mu(:,:,z,:), f(:,:,z,:), c(:,:,z,:));
                end
            end
        end
        
    end
end

function ll = onMemoryProba(mu, f, c)
    mu    = single(numeric(mu));
    f     = single(numeric(f));
    c     = single(numeric(c));
    msk   = isfinite(f) & isfinite(mu) & c > 0;
    mu    = max(eps('single'), min(1 - eps('single'), mu));
    ll    = sumall(f(msk) .* log(mu(msk)) + (1 - f(msk)) .* log(1 - mu(msk)));
end

function ll = onMemoryLog(a, f, c)
    a     = single(numeric(a));
    f     = single(numeric(f));
    c     = single(numeric(c));
    msk   = isfinite(f) & isfinite(a) & c > 0;
    ll    = sumall(f(msk) .* a(msk) - log(1 + exp(a(msk))));
end