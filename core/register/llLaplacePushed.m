function ll = llLaplacePushed(mu, f, c, varargin)
%__________________________________________________________________________
%
% Log-likelihood of the Laplace matching term
%
% -------------------------------------------------------------------------
%
% FORMAT ll = llLaplacePushed(mu, f, c, (b), ...)
%
% REQUIRED
% --------
% mu   - Template
% f    - Observed image pushed in the template space.
% c    - Pushed voxel count.
% 
% OPTIONAL
% --------
% b    - Noise variance (one per modality).
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
% ll   - Log-likelihood of the Laplace matching term
%__________________________________________________________________________

    % --- Parse inputs
    p = inputParser;
    p.FunctionName = 'llLaplacePushed';
    p.addRequired('mu',  @checkarray);
    p.addRequired('f',   @checkarray);
    p.addRequired('c',   @checkarray);
    p.addOptional('b',   1,  @(X) isscalar(X) || checkarray(X));
    p.addParameter('bb',     struct, @isstruct);
    p.addParameter('loop',   '',    @ischar);
    p.addParameter('par',    false, @isscalar);
    p.addParameter('debug',  false, @isscalar);
    p.parse(mu, f, c, varargin{:});
    b    = p.Results.b;
    bb   = p.Results.bb;
    par  = p.Results.par;
    loop = p.Results.loop;
    
    if p.Results.debug, fprintf('* llLaplacePushed\n'); end;
    
    % --- Optimise parallelisation and splitting schemes
    [par, loop] = autoParLoop(par, loop, isa(mu, 'file_array'), size(mu, 3), size(mu, 4));
    
    % --- Prepare sigma
    nc = size(mu, 4);
    if length(b) == 1
        b = b * ones(1, nc);
    end
    
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
    dim  = [numel(bb.x)  numel(bb.y) numel(bb.z) nc];
    dlat = dim(1:3);
    
    % --- Initialise
    ll = 0;
    
    % --- No loop
    if strcmpi(loop, 'none')
        if p.Results.debug, fprintf('   - No loop\n'); end;
        ll = ll + onMemory(mu(bb.x,bb.y,bb.z,:), f, c, s);
    
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
                ll = ll + onMemory(mu(bb.x,bb.y,bb.z,k), f(:,:,:,k), c, s(k));
            end
        elseif isa(mu, 'file_array') && isa(f, 'file_array')
            parfor (k=1:nc, par)
                ll = ll + onMemory(slicevol(mu, {bb.x, bb.y, bb.z, k}), slicevol(f, k, 4), c, s(k));
            end
        elseif isa(mu, 'file_array')
            parfor (k=1:nc, par)
                ll = ll + onMemory(slicevol(mu, {bb.x, bb.y, bb.z, k}), f(:,:,:,k), c, s(k));
            end
        elseif isa(f, 'file_array')
            mu = mu(bb.x, bb.y, bb.z,:);
            parfor (k=1:nc, par)
                ll = ll + onMemory(mu(:,:,:,k), slicevol(f, k, 4), c, s(k));
            end
        else
            mu = mu(bb.x, bb.y, bb.z,:);
            parfor (k=1:nc, par)
                ll = ll + onMemory(mu(:,:,:,k), f(:,:,:,k), c, s(k));
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
                ll = ll + onMemory(mu(bb.x,bb.y,bb.z(1)+z-1,:), f(:,:,z,:), c(:,:,z,:), s);
            end
        elseif isa(mu, 'file_array') && isa(f, 'file_array')
            parfor (z=1:dlat(3), par)
                ll = ll + onMemory(slicevol(mu, {bb.x, bb.y, bb.z(1)+z-1}), slicevol(f, z, 3), slicevol(c, z, 3), b);
            end
        elseif isa(mu, 'file_array')
            parfor (z=1:dlat(3), par)
                ll = ll + onMemory(slicevol(mu, {bb.x, bb.y, bb.z(1)+z-1}), f(:,:,z,:), c(:,:,z,:), b);
            end
        elseif isa(f, 'file_array')
            mu = mu(bb.x, bb.y, bb.z,:);
            parfor (z=1:dlat(3), par)
                ll = ll + onMemory(mu(:,:,z,:), slicevol(f, z, 3), slicevol(c, z, 3), b);
            end
        else
            mu = mu(bb.x, bb.y, bb.z,:);
            parfor (z=1:dlat(3), par)
                ll = ll + onMemory(mu(:,:,z,:), f(:,:,z,:), c(:,:,z,:), b);
            end
        end
        
    end
end

function ll = onMemory(mu, f, c, b)
    nc    = size(mu, 4);
    mu    = single(numeric(mu));
    f     = single(numeric(f));
    c     = single(numeric(c));
    msk   = isfinite(f) & isfinite(mu);
    msk   = bsxfun(@and, msk, c > 0);
    count = sum(sum(sum(msk, 3), 2), 1);
    ll    = 0;
    for k=1:nc
        ll = ll - count(k) * log(2 * b(k)) ...
                - 1 / b(k) * sumall((mu(msk(:,:,:,k)) - f(msk(:,:,:,k))).^2);
    end
end
