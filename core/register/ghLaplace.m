function [g, h, htype] = ghLaplace(mu, f, varargin)
%__________________________________________________________________________
%
% Gradient & Hessian of the **negative** log-likelihood of the laplace 
% matching term w.r.t. changes in the initial velocity.
%
%--------------------------------------------------------------------------
%
% FORMAT [g, (h, htype)] = ghLaplace(mu, f, (b), (gmu), ...)
% REQUIRED
% --------
% mu    - Template
% f     - Observed image pushed in the template space.
%
% OPTIONAL
% --------
% b     - Noise parameter (one per modality).
% gmu   - Template spatial gradients.
%
% KEYWORD ARGUMENTS
% -----------------
% count - Pushed voxel count.
% bb    - Bounding box (if different between template and pushed image)
% loop  - Specify how to split data processing
%         ('component', 'slice' or 'none' [default])
% par   - If true, parallelise processing [default: false]
%
% OUTPUT
% --------
% g     - Gradient
% h     - Hessian
% htype - Shape of the hessian ('diagonal', 'symtensor')
%
%--------------------------------------------------------------------------
%
% Let us remind that the matching term is of the form
%   E = sum c.F(f/c, mu)
% Hence,
%   dE = c.F'(f/c, mu) * dmu
%   d2E = tr(dmu) * c.F''(f/c, mu) * dmu
%
% If gmu is provided
%   Returns gradiant/hessian of the matching term w.r.t. changes in v.
%   size(g) = [nx ny nz 3]
%   size(h) = [nx ny nz 6]
%   htype   = 'symtensor'
% Else
%   Returns only gradients and hessian of the function F of f and mu
%   size(g) = [nx ny nz nc]
%   size(h) = [nx ny nz nc]
%   htype   = 'diagonal'
%__________________________________________________________________________

    % --- Parse inputs
    p = inputParser;
    p.FunctionName = 'ghLaplace';
    p.addRequired('mu',  @checkarray);
    p.addRequired('f',   @checkarray);
    p.addRequired('c',   @checkarray);
    p.addOptional('b',   1,  @checkarray);
    p.addOptional('gmu', []);
    p.addParameter('count',  [],     @(X) isnumeric(X) || isa(X, 'file_array'));
    p.addParameter('bb',     struct, @isstruct);
    p.addParameter('loop',   '',    @ischar);
    p.addParameter('par',    false, @isscalar);
    p.addParameter('output', []);
    p.addParameter('debug',  false, @isscalar);
    p.parse(mu, f, varargin{:});
    gmu  = p.Results.gmu;
    c    = p.Results.count;
    bb   = p.Results.bb;
    par  = p.Results.par;
    loop = p.Results.loop;
    
    if p.Results.debug, fprintf('* ghLaplace\n'); end
    
    % --- Optimise parallelisation and splitting schemes
    [par, loop] = autoParLoop(par, loop, isa(mu, 'file_array'), size(mu, 3), size(mu, 4));
    
    % --- Prepare sigma
    nc = size(mu, 4);
    b = p.Results.b;
    if length(b) == 1
        b = b * ones(1, nc);
    end
    c = single(numeric(c));
    
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
    dim = [numel(bb.x)  numel(bb.y) numel(bb.z) size(mu, 4)];
    dlat = dim(1:3);
    nc   = dim(4);
    nvec = size(gmu, 5);
    
    % --- Prepare output
    output = p.Results.output;
    if ~iscell(output)
        output = {output};
    end
    if numel(output) < 2
        output = [output {[]}];
    end
    if isempty(gmu)
        g = prepareOnDisk(output{1}, [dlat nc]);
        g(:) = 0;
    else
        g = prepareOnDisk(output{1}, [dlat nvec]);
        g(:) = 0;
    end
    if nargout > 1
        if isempty(gmu)
            h = prepareOnDisk(output{2}, [dlat nc]);
            h(:) = 0;
        else
            h = prepareOnDisk(output{2}, [dlat nvec*(nvec+1)/2]);
            h(:) = 0;
        end
    end
    
    % --- No loop
    if strcmpi(loop, 'none')
        if p.Results.debug, fprintf('   - No loop\n'); end;
        if isempty(gmu)
            if nargout > 1
                [g(:,:,:,:), h(:,:,:,:)] = onMemory(mu, f, c, b);
            else
                g(:,:,:,:) = onMemory(mu, f, c, b);
            end
        else
            if nargout > 1
                [g(:,:,:,:), h(:,:,:,:)] = onMemory(mu, f, c, b, gmu);
            else
                g(:,:,:,:) = onMemory(mu, f, c, b, gmu);
            end
        end
    
    % --- Loop on components
    elseif strcmpi(loop, 'component')
        % - Parallelised
        if p.Results.debug
            if par > 0
                fprintf('   - Parallelise on components\n'); 
            else
                fprintf('   - Serialise on components\n'); 
            end
        end
        if isempty(gmu)
            if nargout > 1
                if ~par
                    for k=1:nc
                        [g(:,:,:,k), h(:,:,:,k)] = onMemory(mu(bb.x,bb.y,bb.z,k), f(:,:,:,k), c, b(k));
                    end
                elseif isa(mu, 'file_array') && isa(f, 'file_array')
                    parfor (k=1:nc, par)
                        [g(:,:,:,k), h(:,:,:,k)] = onMemory(slicevol(mu, {bb.x, bb.y, bb.z, k}), slicevol(f, k, 4), c, b(k));
                    end
                elseif isa(mu, 'file_array')
                    parfor (k=1:nc, par)
                        [g(:,:,:,k), h(:,:,:,k)] = onMemory(slicevol(mu, {bb.x, bb.y, bb.z, k}), f(:,:,:,k), c, b(k));
                    end
                elseif isa(f, 'file_array')
                    mu = mu(bb.x,bb.y,bb.z,:);
                    parfor (k=1:nc, par)
                        [g(:,:,:,k), h(:,:,:,k)] = onMemory(mu(:,:,:,k), slicevol(f, k, 4), c, b(k));
                    end
                else
                    mu = mu(bb.x,bb.y,bb.z,:);
                    parfor (k=1:nc, par)
                        [g(:,:,:,k), h(:,:,:,k)] = onMemory(mu(:,:,:,k), f(:,:,:,k), c, b(k));
                    end
                end
            else
                if ~par
                    for k=1:nc
                        g(:,:,:,k) = onMemory(mu(bb.x,bb.y,bb.z,k), f(:,:,:,k), c, b(k));
                    end
                elseif isa(mu, 'file_array') && isa(f, 'file_array')
                    parfor (k=1:nc, par)
                        g(:,:,:,k) = onMemory(slicevol(mu, {bb.x, bb.y, bb.z, k}), slicevol(f, k, 4), c, b(k));
                    end
                elseif isa(mu, 'file_array')
                    parfor (k=1:nc, par)
                        g(:,:,:,k) = onMemory(slicevol(mu, {bb.x, bb.y, bb.z, k}), f(:,:,:,k), c, b(k));
                    end
                elseif isa(f, 'file_array')
                    mu = mu(bb.x,bb.y,bb.z,:);
                    parfor (k=1:nc, par)
                        g(:,:,:,k) = onMemory(mu(:,:,:,k), slicevol(f, k, 4), c, b(k));
                    end
                else
                    mu = mu(bb.x,bb.y,bb.z,:);
                    parfor (k=1:nc, par)
                        g(:,:,:,k) = onMemory(mu(:,:,:,k), f(:,:,:,k), c, b(k));
                    end
                end
            end
        else
            if nargout > 1
                g = numeric(g);
                h = numeric(h);
                if ~par
                    for k=1:nc
                        [g1, h1] = onMemory(mu(bb.x,bb.y,bb.z,k), f(:,:,:,k), c, b(:,k), gmu(bb.x,bb.y,bb.z,k,:));
                        g = g + g1;
                        h = h + h1;
                    end
                elseif isa(mu, 'file_array') && isa(gmu, 'file_array') && isa(f, 'file_array')
                    parfor (k=1:nc, par)
                        [g1, h1] = onMemory(slicevol(mu, {bb.x, bb.y, bb.z, k}), slicevol(f, k, 4), c, b(:,k), slicevol(gmu, {bb.x, bb.y, bb.z, k}));
                        g = g + g1;
                        h = h + h1;
                    end
                elseif isa(mu, 'file_array') && isa(gmu, 'file_array')
                    parfor (k=1:nc, par)
                        [g1, h1] = onMemory(slicevol(mu, {bb.x, bb.y, bb.z, k}), f(:,:,:,k), c, b(:,k), slicevol(gmu, {bb.x, bb.y, bb.z, k}));
                        g = g + g1;
                        h = h + h1;
                    end
                elseif isa(mu, 'file_array') && isa(f, 'file_array')
                    gmu = mu(bb.x,bb.y,bb.z,:,:);
                    parfor (k=1:nc, par)
                        [g1, h1] = onMemory(slicevol(mu, {bb.x, bb.y, bb.z, k}), slicevol(f, k, 4), c, b(:,k), gmu(:,:,:,k,:));
                        g = g + g1;
                        h = h + h1;
                    end
                elseif isa(gmu, 'file_array') && isa(f, 'file_array')
                    mu = mu(bb.x,bb.y,bb.z,:);
                    parfor (k=1:nc, par)
                        [g1, h1] = onMemory(mu(:,:,:,k), slicevol(f, k, 4), c, b(:,k), slicevol(gmu, {bb.x, bb.y, bb.z, k}));
                        g = g + g1;
                        h = h + h1;
                    end
                elseif isa(mu, 'file_array')
                    gmu = mu(bb.x,bb.y,bb.z,:,:);
                    parfor (k=1:nc, par)
                        [g1, h1] = onMemory(slicevol(mu, {bb.x, bb.y, bb.z, k}), f(:,:,:,k), c, b(:,k), gmu(:,:,:,k,:));
                        g = g + g1;
                        h = h + h1;
                    end
                elseif isa(gmu, 'file_array')
                    mu = mu(bb.x,bb.y,bb.z,:);
                    parfor (k=1:nc, par)
                        [g1, h1] = onMemory(mu(:,:,:,k), f(:,:,:,k), c, b(:,k), slicevol(gmu, {bb.x, bb.y, bb.z, k}));
                        g = g + g1;
                        h = h + h1;
                    end
                elseif isa(f, 'file_array')
                    mu = mu(bb.x,bb.y,bb.z,:);
                    gmu = mu(bb.x,bb.y,bb.z,:,:);
                    parfor (k=1:nc, par)
                        [g1, h1] = onMemory(mu(:,:,:,k), slicevol(f, k, 4), c, b(:,k), gmu(:,:,:,k,:));
                        g = g + g1;
                        h = h + h1;
                    end
                else
                    mu  = mu(bb.x,bb.y,bb.z,:);
                    gmu = mu(bb.x,bb.y,bb.z,:,:);
                    parfor (k=1:nc, par)
                        [g1, h1] = onMemory(mu(:,:,:,k), f(:,:,:,k), c, b(:,k), gmu(:,:,:,k,:));
                        g = g + g1;
                        h = h + h1;
                    end
                end
            else
                g = numeric(g);
                if ~par
                    for k=1:nc
                        g1 = onMemory(mu(bb.x,bb.y,bb.z,k), f(:,:,:,k), c, b(:,k), gmu(bb.x,bb.y,bb.z,k,:));
                        g = g + g1;
                    end
                elseif isa(mu, 'file_array') && isa(gmu, 'file_array') && isa(f, 'file_array')
                    parfor (k=1:nc, par)
                        g1 = onMemory(slicevol(mu, {bb.x, bb.y, bb.z, k}), slicevol(f, k, 4), c, b(:,k), slicevol(gmu, {bb.x, bb.y, bb.z, k}));
                        g = g + g1;
                    end
                elseif isa(mu, 'file_array') && isa(gmu, 'file_array')
                    parfor (k=1:nc, par)
                        g1 = onMemory(slicevol(mu, {bb.x, bb.y, bb.z, k}), f(:,:,:,k), c, b(:,k), slicevol(gmu, {bb.x, bb.y, bb.z, k}));
                        g = g + g1;
                    end
                elseif isa(mu, 'file_array') && isa(f, 'file_array')
                    gmu = mu(bb.x,bb.y,bb.z,:,:);
                    parfor (k=1:nc, par)
                        g1 = onMemory(slicevol(mu, {bb.x, bb.y, bb.z, k}), slicevol(f, k, 4), c, b(:,k), gmu(:,:,:,k,:));
                        g = g + g1;
                    end
                elseif isa(gmu, 'file_array') && isa(f, 'file_array')
                    mu  = mu(bb.x,bb.y,bb.z,:);
                    parfor (k=1:nc, par)
                        g1 = onMemory(mu(:,:,:,k), slicevol(f, k, 4), c, b(:,k), slicevol(gmu, {bb.x, bb.y, bb.z, k}));
                        g = g + g1;
                    end
                elseif isa(mu, 'file_array')
                    gmu = mu(bb.x,bb.y,bb.z,:,:);
                    parfor (k=1:nc, par)
                        g1 = onMemory(slicevol(mu, {bb.x, bb.y, bb.z, k}), f(:,:,:,k), c, b(:,k), gmu(:,:,:,k,:));
                        g = g + g1;
                    end
                elseif isa(gmu, 'file_array')
                    mu  = mu(bb.x,bb.y,bb.z,:);
                    parfor (k=1:nc, par)
                        g1 = onMemory(mu(:,:,:,k), f(:,:,:,k), c, b(:,k), slicevol(gmu, {bb.x, bb.y, bb.z, k}));
                        g = g + g1;
                    end
                elseif isa(f, 'file_array')
                    mu  = mu(bb.x,bb.y,bb.z,:);
                    gmu = mu(bb.x,bb.y,bb.z,:,:);
                    parfor (k=1:nc, par)
                        g1 = onMemory(mu(:,:,:,k), slicevol(f, k, 4), c, b(:,k), gmu(:,:,:,k,:));
                        g = g + g1;
                    end
                else
                    mu  = mu(bb.x,bb.y,bb.z,:);
                    gmu = mu(bb.x,bb.y,bb.z,:,:);
                    parfor (k=1:nc, par)
                        g1 = onMemory(mu(:,:,:,k), f(:,:,:,k), c, b(:,k), gmu(:,:,:,k,:));
                        g = g + g1;
                    end
                end
            end
        end
        
    % --- Loop on slices
    elseif strcmpi(loop, 'slice')
        % - Parallelised
        if p.Results.debug
            if par > 0
                fprintf('   - Parallelise on slices\n'); 
            else
                fprintf('   - Serialise on slices\n'); 
            end
        end
        if isempty(gmu)
            if nargout > 1
                if ~par
                    for z=1:dlat(3)
                        [g1, h1] = onMemory(mu(bb.x,bb.y,bb.z(1)+z-1,:), f(:,:,z,:), c(:,:,z), b);
                        g(:,:,z,:) = g(:,:,z,:) + g1;
                        h(:,:,z,:) = h(:,:,z,:) + h1;
                    end
                elseif isa(mu, 'file_array') && isa(f, 'file_array')
                    parfor (z=1:dlat(3), par)
                        [g1, h1] = onMemory(slicevol(mu, {bb.x,bb.y,bb.z(1)+z-1}), slicevol(f, z, 3), slicevol(c, z, 3), b);
                        g(:,:,z,:) = g(:,:,z,:) + g1;
                        h(:,:,z,:) = h(:,:,z,:) + h1;
                    end
                elseif isa(mu, 'file_array')
                    parfor (z=1:dlat(3), par)
                        [g1, h1] = onMemory(slicevol(mu, {bb.x,bb.y,bb.z(1)+z-1}), f(:,:,z,:), c(:,:,z), b);
                        g(:,:,z,:) = g(:,:,z,:) + g1;
                        h(:,:,z,:) = h(:,:,z,:) + h1;
                    end
                elseif isa(f, 'file_array')
                    mu = mu(bb.x, bb.y, bb.z, :);
                    parfor (z=1:dlat(3), par)
                        [g1, h1] = onMemory(mu(:,:,z,:), slicevol(f, z, 3), slicevol(c, z, 3), b);
                        g(:,:,z,:) = g(:,:,z,:) + g1;
                        h(:,:,z,:) = h(:,:,z,:) + h1;
                    end
                else
                    mu = mu(bb.x, bb.y, bb.z, :);
                    parfor (z=1:dlat(3), par)
                        [g1, h1] = onMemory(mu(:,:,z,:), f(:,:,z,:), c(:,:,z), b);
                        g(:,:,z,:) = g(:,:,z,:) + g1;
                        h(:,:,z,:) = h(:,:,z,:) + h1;
                    end
                end
             else
                if ~par
                    for z=1:dlat(3)
                        g1 = onMemory(mu(bb.x,bb.y,bb.z(1)+z-1,:), f(:,:,z,:), c(:,:,z), b);
                        g(:,:,z,:) = g(:,:,z,:) + g1;
                    end
                elseif isa(mu, 'file_array') && isa(f, 'file_array')
                    parfor (z=1:dlat(3), par)
                        g1 = onMemory(slicevol(mu, {bb.x,bb.y,bb.z(1)+z-1}), slicevol(f, z, 3), slicevol(c, z, 3), b);
                        g(:,:,z,:) = g(:,:,z,:) + g1;
                    end
                elseif isa(mu, 'file_array')
                    parfor (z=1:dlat(3), par)
                        g1 = onMemory(slicevol(mu, {bb.x,bb.y,bb.z(1)+z-1}), f(:,:,z,:), c(:,:,z), b);
                        g(:,:,z,:) = g(:,:,z,:) + g1;
                    end
                elseif isa(f, 'file_array')
                    mu = mu(bb.x, bb.y, bb.z, :);
                    parfor (z=1:dlat(3), par)
                        g1 = onMemory(mu(:,:,z,:), slicevol(f, z, 3), slicevol(c, z, 3), b);
                        g(:,:,z,:) = g(:,:,z,:) + g1;
                    end
                else
                    mu = mu(bb.x, bb.y, bb.z, :);
                    parfor (z=1:dlat(3), par)
                        g1 = onMemory(mu(:,:,z,:), f(:,:,z,:), c(:,:,z), b);
                        g(:,:,z,:) = g(:,:,z,:) + g1;
                    end
                end
            end
        else
            if nargout > 1
                g = numeric(g);
                h = numeric(h);
                if ~par
                    for z=1:dlat(3)
                        [g1, h1] = onMemory(mu(bb.x,bb.y,bb.z(1)+z-1,:), f(:,:,z,:), c(:,:,z), b, gmu(bb.x,bb.y,bb.z(1)+z-1,:,:));
                        g(:,:,z,:) = g(:,:,z,:) + g1;
                        h(:,:,z,:) = h(:,:,z,:) + h1;
                    end
                elseif isa(mu, 'file_array') && isa(gmu, 'file_array') && isa(f, 'file_array')
                    parfor (z=1:dlat(3), par)
                        [g1, h1] = onMemory(slicevol(mu, {bb.x,bb.y,bb.z(1)+z-1}), slicevol(f, z, 3), slicevol(c, z, 3), b, slicevol(gmu, {bb.x,bb.y,bb.z(1)+z-1}));
                        g(:,:,z,:) = g(:,:,z,:) + g1;
                        h(:,:,z,:) = h(:,:,z,:) + h1;
                    end
                elseif isa(mu, 'file_array') && isa(gmu, 'file_array')
                    parfor (z=1:dlat(3), par)
                        [g1, h1] = onMemory(slicevol(mu, {bb.x,bb.y,bb.z(1)+z-1}), f(:,:,z,:), c(:,:,z), b, slicevol(gmu, {bb.x,bb.y,bb.z(1)+z-1}));
                        g(:,:,z,:) = g(:,:,z,:) + g1;
                        h(:,:,z,:) = h(:,:,z,:) + h1;
                    end
                elseif isa(mu, 'file_array') && isa(f, 'file_array')
                    gmu = gmu(bb.x,bb.y,bb.z,:,:);
                    parfor (z=1:dlat(3), par)
                        [g1, h1] = onMemory(slicevol(mu, {bb.x,bb.y,bb.z(1)+z-1}), slicevol(f, z, 3), slicevol(c, z, 3), b, gmu(:,:,z,:,:));
                        g(:,:,z,:) = g(:,:,z,:) + g1;
                        h(:,:,z,:) = h(:,:,z,:) + h1;
                    end
                elseif isa(gmu, 'file_array') && isa(f, 'file_array')
                    mu  = mu(bb.x,bb.y,bb.z,:);
                    gmu = gmu(bb.x,bb.y,bb.z,:,:);
                    parfor (z=1:dlat(3), par)
                        [g1, h1] = onMemory(mu(:,:,z,:), slicevol(f, z, 3), slicevol(c, z, 3), b, slicevol(gmu, {bb.x,bb.y,bb.z(1)+z-1}));
                        g(:,:,z,:) = g(:,:,z,:) + g1;
                        h(:,:,z,:) = h(:,:,z,:) + h1;
                    end
                elseif isa(mu, 'file_array')
                    gmu = gmu(bb.x,bb.y,bb.z,:,:);
                    parfor (z=1:dlat(3), par)
                        [g1, h1] = onMemory(slicevol(mu, {bb.x,bb.y,bb.z(1)+z-1}), f(:,:,z,:), c(:,:,z), b, gmu(:,:,z,:,:));
                        g(:,:,z,:) = g(:,:,z,:) + g1;
                        h(:,:,z,:) = h(:,:,z,:) + h1;
                    end
                elseif isa(gmu, 'file_array')
                    mu  = mu(bb.x,bb.y,bb.z,:);
                    parfor (z=1:dlat(3), par)
                        [g1, h1] = onMemory(mu(:,:,z,:), f(:,:,z,:), c(:,:,z), b, slicevol(gmu, {bb.x,bb.y,bb.z(1)+z-1}));
                        g(:,:,z,:) = g(:,:,z,:) + g1;
                        h(:,:,z,:) = h(:,:,z,:) + h1;
                    end
                elseif isa(f, 'file_array')
                    mu  = mu(bb.x,bb.y,bb.z,:);
                    gmu = gmu(bb.x,bb.y,bb.z,:,:);
                    parfor (z=1:dlat(3), par)
                        [g1, h1] = onMemory(mu(:,:,z,:), slicevol(f, z, 3), slicevol(c, z, 3), b, gmu(:,:,z,:,:));
                        g(:,:,z,:) = g(:,:,z,:) + g1;
                        h(:,:,z,:) = h(:,:,z,:) + h1;
                    end
                else
                    mu  = mu(bb.x,bb.y,bb.z,:);
                    gmu = gmu(bb.x,bb.y,bb.z,:,:);
                    parfor (z=1:dlat(3), par)
                        [g1, h1] = onMemory(mu(:,:,z,:), f(:,:,z,:), c(:,:,z), b, gmu(:,:,z,:,:));
                        g(:,:,z,:) = g(:,:,z,:) + g1;
                        h(:,:,z,:) = h(:,:,z,:) + h1;
                    end
                end
            else
                g = numeric(g);
                if ~par
                    for z=1:dlat(3)
                        g1 = onMemory(mu(bb.x,bb.y,bb.z(1)+z-1,:), f(:,:,z,:), c(:,:,z), b, gmu(bb.x,bb.y,bb.z(1)+z-1,:,:));
                        g(:,:,z,:) = g(:,:,z,:) + g1;
                    end
                elseif isa(mu, 'file_array') && isa(gmu, 'file_array') && isa(f, 'file_array')
                    parfor (z=1:dlat(3), par)
                        g1 = onMemory(slicevol(mu, {bb.x,bb.y,bb.z(1)+z-1}), slicevol(f, z, 3), slicevol(c, z, 3), b, slicevol(gmu, {bb.x,bb.y,bb.z(1)+z-1}));
                        g(:,:,z,:) = g(:,:,z,:) + g1;
                    end
                elseif isa(mu, 'file_array') && isa(gmu, 'file_array')
                    parfor (z=1:dlat(3), par)
                        g1 = onMemory(slicevol(mu, {bb.x,bb.y,bb.z(1)+z-1}), f(:,:,z,:), c(:,:,z), b, slicevol(gmu, {bb.x,bb.y,bb.z(1)+z-1}));
                        g(:,:,z,:) = g(:,:,z,:) + g1;
                    end
                elseif isa(mu, 'file_array') && isa(f, 'file_array')
                    gmu = gmu(bb.x,bb.y,bb.z,:,:);
                    parfor (z=1:dlat(3), par)
                        g1 = onMemory(slicevol(mu, {bb.x,bb.y,bb.z(1)+z-1}), slicevol(f, z, 3), slicevol(c, z, 3), b, gmu(:,:,z,:,:));
                        g(:,:,z,:) = g(:,:,z,:) + g1;
                    end
                elseif isa(gmu, 'file_array') && isa(f, 'file_array')
                    mu  = mu(bb.x,bb.y,bb.z,:);
                    parfor (z=1:dlat(3), par)
                        g1 = onMemory(mu(:,:,z,:), slicevol(f, z, 3), slicevol(c, z, 3), b, slicevol(gmu, {bb.x,bb.y,bb.z(1)+z-1}));
                        g(:,:,z,:) = g(:,:,z,:) + g1;
                    end
                elseif isa(mu, 'file_array')
                    gmu = gmu(bb.x,bb.y,bb.z,:,:);
                    parfor (z=1:dlat(3), par)
                        g1 = onMemory(slicevol(mu, {bb.x,bb.y,bb.z(1)+z-1}), f(:,:,z,:), c(:,:,z), b, gmu(:,:,z,:,:));
                        g(:,:,z,:) = g(:,:,z,:) + g1;
                    end
                elseif isa(gmu, 'file_array')
                    mu  = mu(bb.x,bb.y,bb.z,:);
                    parfor (z=1:dlat(3), par)
                        g1 = onMemory(mu(:,:,z,:), f(:,:,z,:), c(:,:,z), b, slicevol(gmu, {bb.x,bb.y,bb.z(1)+z-1}));
                        g(:,:,z,:) = g(:,:,z,:) + g1;
                    end
                elseif isa(f, 'file_array')
                    mu  = mu(bb.x,bb.y,bb.z,:);
                    gmu = gmu(bb.x,bb.y,bb.z,:,:);
                    parfor (z=1:dlat(3), par)
                        g1 = onMemory(mu(:,:,z,:), slicevol(f, z, 3), slicevol(c, z, 3), b, gmu(:,:,z,:,:));
                        g(:,:,z,:) = g(:,:,z,:) + g1;
                    end
                else
                    mu  = mu(bb.x,bb.y,bb.z,:);
                    gmu = gmu(bb.x,bb.y,bb.z,:,:);
                    parfor (z=1:dlat(3), par)
                        g1 = onMemory(mu(:,:,z,:), f(:,:,z,:), c(:,:,z), b, gmu(:,:,z,:,:));
                        g(:,:,z,:) = g(:,:,z,:) + g1;
                    end
                end
            end
        end
        
    end
    
    % --- Set hessian type
    if isempty(gmu)
        htype = 'diagonal';
    else
        htype = 'symtensor';
    end
    
    % --- Write on disk
    if hessian
        if ~isempty(output{1})
            h = saveOnDisk(output{1}, h, 'name', 'h');
        end
        g = [];
        [g, h, htype] = deal(h, htype, g);
    elseif nargin > 1
        if ~isempty(output{1})
            g = saveOnDisk(output{1}, g, 'name', 'g');
        end
        if ~isempty(output{2})
            h = saveOnDisk(output{2}, h, 'name', 'h');
        end
    else
        if ~isempty(output{1})
            g = saveOnDisk(output{1}, g, 'name', 'g');
        end
    end
    
end

% Actual gradient and hessian computation
function [g, h] = onMemory(mu, f, c, b, gmu)
    
    if nargin < 5
        gmu = [];
    end
    lat = [size(mu) 1 1];
    lat = lat(1:3);
    nc  = size(mu, 4);
    
    mu  = single(numeric(mu));
    f   = single(numeric(f));
    c   = single(numeric(c));
    if isempty(c)
        c = single(1);
    end
    b   = reshape(1./b, [1 1 1 nc]);
    
    if ~hessian
        g  = sign(bsxfun(@times, c, mu) - f);
        g  = bsxfun(@times, b, g);
        if ~isempty(gmu)
            g = -spm_matcomp('Pointwise', gmu, g, 't');
        end
        g(~isfinite(g)) = 0;
    else
        g= [];
    end
    if nargout > 1 || hessian
        h = bsxfun(@times, b, c);
        if ~isempty(gmu)
            nvec = size(gmu, 5);
            [ind, k] = spm_matcomp('SymIndices', nvec, 'n');
            hh = h;
            h = zeros([lat k], 'like', h);
            for d=1:nvec
                for l=d:nvec
                    h(:,:,:,ind(d,l)) = h(:,:,:,ind(d,l)) ...
                        + sum(hh .* gmu(:,:,:,:,d) .* gmu(:,:,:,:,l), 4);
                end
            end
        end
        h(~isfinite(h)) = 0;
    end
end