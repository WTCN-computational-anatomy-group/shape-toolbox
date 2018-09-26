function [g, h, htype] = ghNormal(mu, f, varargin)
%__________________________________________________________________________
%
% Gradient & Hessian of the **negative** log-likelihood of the normal 
% matching term w.r.t. changes in template value or initial velocity.
%
%--------------------------------------------------------------------------
%
% FORMAT [(g), (h, htype)] = ghNormal(mu, f, (s), (ga), ...)
%
% REQUIRED
% --------
% mu    - Reconstructed probability template. [mx my mz k]
% f     - Observed image. [mx my mz k]
%
% Note: we need either mu warped to image space or f pushed to template
%       space.
% 
% OPTIONAL
% --------
% s     - Noise variance (one per modality) [1]
% ga    - Spatial gradients of the template coefficients.
%         > If ga provided, but not ipsi, compute pointwise gM * gA
%           In this case, ga size is [mx my mz k 3]
%
% KEYWORD ARGUMENTS
% -----------------
% ipsi    - Inverse (subj to template) warp [mx my mz 3]
%           > If provided on top of ga, compute push(gM) * gA
%             In this case, ga size is [nx ny nz k 3]
% lat     - Output lattice size (not needed if ga provided)
% circ    - (Push only) Boundary conditions for pushing [1]
% count   - Pushed voxel count.
% bb      - Bounding box (if different between template and pushed image)
% hessian - Compute only hessian (not gradient) [false]
% loop    - Specify how to split data processing
%           ('component', 'slice' or 'none') ['none']
% par     - If true, parallelise processing  [false]
%
% OUTPUT
% ------
% g     - Gradient
%         > If ga and ipsi: [nx ny nz 3]
%         > Else if ga:     [mx my mz 3]
%         > Else if ipsi:   [nx ny nz k]
%         > Else:           [mx my mz k]
% h     - Hessian
%         > If ga and ipsi: [nx ny nz 6]
%         > Else if ga:     [mx my mz 6]
%         > Else if ipsi:   [nx ny nz k]
%         > Else:           [mx my mz k]
% htype - Shape of the hessian
%         > If ga:          'symtensor'
%         > Else:           'diagonal'
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
    p.FunctionName = 'ghNormal';
    p.addRequired('mu',  @checkarray);
    p.addRequired('f',   @checkarray);
    p.addOptional('s',   1,  @checkarray);
    p.addOptional('gmu', []);
    p.addParameter('ipsi',   []);
    p.addParameter('lat',    []);
    p.addParameter('circ',    1);
    p.addParameter('count',  [],     @(X) isnumeric(X) || isa(X, 'file_array'));
    p.addParameter('bb',     struct, @isstruct);
    p.addParameter('hessian', false, @islogical);
    p.addParameter('loop',   '',    @ischar);
    p.addParameter('par',    false, @isscalar);
    p.addParameter('output', []);
    p.addParameter('debug',  false, @isscalar);
    p.parse(mu, f, varargin{:});
    gmu  = p.Results.gmu;
    ipsi = p.Results.ipsi;
    lat  = p.Results.lat;
    c    = p.Results.count;
    bb   = p.Results.bb;
    hessian = p.Results.hessian;
    par  = p.Results.par;
    loop = p.Results.loop;
    
    if p.Results.debug, fprintf('* ghNormal\n'); end
    
    % --- Optimise parallelisation and splitting schemes
    [par, loop] = autoParLoop(par, loop, isa(mu, 'file_array'), size(mu, 3), size(mu, 4));
    
    % --- Prepare sigma
    nc = size(mu, 4);
    s = p.Results.s;
    if length(s) == 1
        s = s * ones(1, nc);
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
    bbx = bb.x;
    bby = bb.y;
    bbz = bb.z;
    oz = bbz(1) - 1;
    
    % --- Read dimensions
    dim = [numel(bbx)  numel(bby) numel(bbz) size(mu, 4)];
    dlat = dim(1:3);
    nc   = dim(4);
    nvec = size(gmu, 5);
    
    % --- Switch between cases [push(g * pull(gmu))] and [push(g) * gmu]
    if ~isempty(ipsi)
        % We will used gmu for "gradient in image space"
        %         and  ga  for "gradient in template space"
        [ga, gmu] = deal(gmu, []);
        if isempty(lat)
            if ~isempty(ga)
                lat = [size(ga) 1];
            else
                lat = dlat;
            end
        end
        lat = lat(1:3);
    else
        ga = [];
    end
    
    % --- Prepare output
    output = p.Results.output;
    if ~iscell(output)
        output = {output};
    end
    if numel(output) < 2
        output = [output {[]}];
    end
    if ~hessian
        if isempty(gmu)
            g    = prepareOnDisk(output{1}, [dlat nc]);
            g(:) = 0;
        else
            g    = prepareOnDisk(output{1}, [dlat nvec]);
            g(:) = 0;
        end
    end
    if hessian || nargout > 1
        if isempty(gmu)
            h    = prepareOnDisk(output{2}, [dlat nc]);
            h(:) = 0;
        else
            h    = prepareOnDisk(output{2}, [dlat nvec*(nvec+1)/2]);
            h(:) = 0;
        end
    end
    
    
    % --- No loop
    if strcmpi(loop, 'none')
        if p.Results.debug, fprintf('   - No loop\n'); end
        if isempty(gmu)
            if hessian
                h(:,:,:,:) = onMemory(mu, f, c, s, [], true);
            elseif nargout > 1
                [g(:,:,:,:), h(:,:,:,:)] = onMemory(mu, f, c, s);
            else
                g(:,:,:,:) = onMemory(mu, f, c, s);
            end
        else
            if hessian
                h(:,:,:,:) = onMemory(mu, f, c, s, gmu, true);
            elseif nargout > 1
                [g(:,:,:,:), h(:,:,:,:)] = onMemory(mu, f, c, s, gmu);
            else
                g(:,:,:,:) = onMemory(mu, f, c, s, gmu);
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
            if hessian
                if ~par
                    for k=1:nc
                        h(:,:,:,k) = onMemory(mu(bbx,bby,bbz,k), f(:,:,:,k), c, s(k), [], true);
                    end
                elseif isa(mu, 'file_array') && isa(f, 'file_array')
                    parfor (k=1:nc, par)
                        h(:,:,:,k) = onMemory(slicevol(mu, {bbx, bby, bbz, k}), slicevol(f, k, 4), c, s(k), [], true);
                    end
                elseif isa(mu, 'file_array')
                    parfor (k=1:nc, par)
                        h(:,:,:,k) = onMemory(slicevol(mu, {bbx, bby, bbz, k}), f(:,:,:,k), c, s(k), [], true);
                    end
                elseif isa(f, 'file_array')
                    mu = mu(bbx,bby,bbz,:);
                    parfor (k=1:nc, par)
                        h(:,:,:,k) = onMemory(mu(:,:,:,k), slicevol(f, k, 4), c, s(k), [], true);
                    end
                else
                    mu = mu(bbx,bby,bbz,:);
                    parfor (k=1:nc, par)
                        h(:,:,:,k) = onMemory(mu(:,:,:,k), f(:,:,:,k), c, s(k), [], true);
                    end
                end
            elseif nargout > 1
                if ~par
                    for k=1:nc
                        [g(:,:,:,k), h(:,:,:,k)] = onMemory(mu(bbx,bby,bbz,k), f(:,:,:,k), c, s(k));
                    end
                elseif isa(mu, 'file_array') && isa(f, 'file_array')
                    parfor (k=1:nc, par)
                        [g(:,:,:,k), h(:,:,:,k)] = onMemory(slicevol(mu, {bbx, bby, bbz, k}), slicevol(f, k, 4), c, s(k));
                    end
                elseif isa(mu, 'file_array')
                    parfor (k=1:nc, par)
                        [g(:,:,:,k), h(:,:,:,k)] = onMemory(slicevol(mu, {bbx, bby, bbz, k}), f(:,:,:,k), c, s(k));
                    end
                elseif isa(f, 'file_array')
                    mu = mu(bbx,bby,bbz,:);
                    parfor (k=1:nc, par)
                        [g(:,:,:,k), h(:,:,:,k)] = onMemory(mu(:,:,:,k), slicevol(f, k, 4), c, s(k));
                    end
                else
                    mu = mu(bbx,bby,bbz,:);
                    parfor (k=1:nc, par)
                        [g(:,:,:,k), h(:,:,:,k)] = onMemory(mu(:,:,:,k), f(:,:,:,k), c, s(k));
                    end
                end
            else
                if ~par
                    for k=1:nc
                        g(:,:,:,k) = onMemory(mu(bbx,bby,bbz,k), f(:,:,:,k), c, s(k));
                    end
                elseif isa(mu, 'file_array') && isa(f, 'file_array')
                    parfor (k=1:nc, par)
                        g(:,:,:,k) = onMemory(slicevol(mu, {bbx, bby, bbz, k}), slicevol(f, k, 4), c, s(k));
                    end
                elseif isa(mu, 'file_array')
                    parfor (k=1:nc, par)
                        g(:,:,:,k) = onMemory(slicevol(mu, {bbx, bby, bbz, k}), f(:,:,:,k), c, s(k));
                    end
                elseif isa(f, 'file_array')
                    mu = mu(bbx,bby,bbz,:);
                    parfor (k=1:nc, par)
                        g(:,:,:,k) = onMemory(mu(:,:,:,k), slicevol(f, k, 4), c, s(k));
                    end
                else
                    mu = mu(bbx,bby,bbz,:);
                    parfor (k=1:nc, par)
                        g(:,:,:,k) = onMemory(mu(:,:,:,k), f(:,:,:,k), c, s(k));
                    end
                end
            end
        else
            if hessian
                h = numeric(h);
                if ~par
                    for k=1:nc
                        h1 = onMemory(mu(bbx,bby,bbz,k), f(:,:,:,k), c, s(:,k), gmu(bbx,bby,bbz,k,:), true);
                        h = h + h1;
                    end
                elseif isa(mu, 'file_array') && isa(gmu, 'file_array') && isa(f, 'file_array')
                    parfor (k=1:nc, par)
                        h1 = onMemory(slicevol(mu, {bbx, bby, bbz, k}), slicevol(f, k, 4), c, s(:,k), slicevol(gmu, {bbx, bby, bbz, k}), true);
                        h = h + h1;
                    end
                elseif isa(mu, 'file_array') && isa(gmu, 'file_array')
                    parfor (k=1:nc, par)
                        h1 = onMemory(slicevol(mu, {bbx, bby, bbz, k}), f(:,:,:,k), c, s(:,k), slicevol(gmu, {bbx, bby, bbz, k}));
                        h = h + h1;
                    end
                elseif isa(mu, 'file_array') && isa(f, 'file_array')
                    gmu = mu(bbx,bby,bbz,:,:);
                    parfor (k=1:nc, par)
                        h1 = onMemory(slicevol(mu, {bbx, bby, bbz, k}), slicevol(f, k, 4), c, s(:,k), gmu(:,:,:,k,:), true);
                        h = h + h1;
                    end
                elseif isa(gmu, 'file_array') && isa(f, 'file_array')
                    mu = mu(bbx,bby,bbz,:);
                    parfor (k=1:nc, par)
                        h1 = onMemory(mu(:,:,:,k), slicevol(f, k, 4), c, s(:,k), slicevol(gmu, {bbx, bby, bbz, k}), true);
                        h = h + h1;
                    end
                elseif isa(mu, 'file_array')
                    gmu = mu(bbx,bby,bbz,:,:);
                    parfor (k=1:nc, par)
                        h1 = onMemory(slicevol(mu, {bbx, bby, bbz, k}), f(:,:,:,k), c, s(:,k), gmu(:,:,:,k,:), true);
                        h = h + h1;
                    end
                elseif isa(gmu, 'file_array')
                    mu = mu(bbx,bby,bbz,:);
                    parfor (k=1:nc, par)
                        h1 = onMemory(mu(:,:,:,k), f(:,:,:,k), c, s(:,k), slicevol(gmu, {bbx, bby, bbz, k}), true);
                        h = h + h1;
                    end
                elseif isa(f, 'file_array')
                    mu = mu(bbx,bby,bbz,:);
                    gmu = mu(bbx,bby,bbz,:,:);
                    parfor (k=1:nc, par)
                        h1 = onMemory(mu(:,:,:,k), slicevol(f, k, 4), c, s(:,k), gmu(:,:,:,k,:), true);
                        h = h + h1;
                    end
                else
                    mu  = mu(bbx,bby,bbz,:);
                    gmu = mu(bbx,bby,bbz,:,:);
                    parfor (k=1:nc, par)
                        h1 = onMemory(mu(:,:,:,k), f(:,:,:,k), c, s(:,k), gmu(:,:,:,k,:), true);
                        h = h + h1;
                    end
                end
            elseif nargout > 1
                g = numeric(g);
                h = numeric(h);
                if ~par
                    for k=1:nc
                        [g1, h1] = onMemory(mu(bbx,bby,bbz,k), f(:,:,:,k), c, s(:,k), gmu(bbx,bby,bbz,k,:));
                        g = g + g1;
                        h = h + h1;
                    end
                elseif isa(mu, 'file_array') && isa(gmu, 'file_array') && isa(f, 'file_array')
                    parfor (k=1:nc, par)
                        [g1, h1] = onMemory(slicevol(mu, {bbx, bby, bbz, k}), slicevol(f, k, 4), c, s(:,k), slicevol(gmu, {bbx, bby, bbz, k}));
                        g = g + g1;
                        h = h + h1;
                    end
                elseif isa(mu, 'file_array') && isa(gmu, 'file_array')
                    parfor (k=1:nc, par)
                        [g1, h1] = onMemory(slicevol(mu, {bbx, bby, bbz, k}), f(:,:,:,k), c, s(:,k), slicevol(gmu, {bbx, bby, bbz, k}));
                        g = g + g1;
                        h = h + h1;
                    end
                elseif isa(mu, 'file_array') && isa(f, 'file_array')
                    gmu = mu(bbx,bby,bbz,:,:);
                    parfor (k=1:nc, par)
                        [g1, h1] = onMemory(slicevol(mu, {bbx, bby, bbz, k}), slicevol(f, k, 4), c, s(:,k), gmu(:,:,:,k,:));
                        g = g + g1;
                        h = h + h1;
                    end
                elseif isa(gmu, 'file_array') && isa(f, 'file_array')
                    mu = mu(bbx,bby,bbz,:);
                    parfor (k=1:nc, par)
                        [g1, h1] = onMemory(mu(:,:,:,k), slicevol(f, k, 4), c, s(:,k), slicevol(gmu, {bbx, bby, bbz, k}));
                        g = g + g1;
                        h = h + h1;
                    end
                elseif isa(mu, 'file_array')
                    gmu = mu(bbx,bby,bbz,:,:);
                    parfor (k=1:nc, par)
                        [g1, h1] = onMemory(slicevol(mu, {bbx, bby, bbz, k}), f(:,:,:,k), c, s(:,k), gmu(:,:,:,k,:));
                        g = g + g1;
                        h = h + h1;
                    end
                elseif isa(gmu, 'file_array')
                    mu = mu(bbx,bby,bbz,:);
                    parfor (k=1:nc, par)
                        [g1, h1] = onMemory(mu(:,:,:,k), f(:,:,:,k), c, s(:,k), slicevol(gmu, {bbx, bby, bbz, k}));
                        g = g + g1;
                        h = h + h1;
                    end
                elseif isa(f, 'file_array')
                    mu = mu(bbx,bby,bbz,:);
                    gmu = mu(bbx,bby,bbz,:,:);
                    parfor (k=1:nc, par)
                        [g1, h1] = onMemory(mu(:,:,:,k), slicevol(f, k, 4), c, s(:,k), gmu(:,:,:,k,:));
                        g = g + g1;
                        h = h + h1;
                    end
                else
                    mu  = mu(bbx,bby,bbz,:);
                    gmu = mu(bbx,bby,bbz,:,:);
                    parfor (k=1:nc, par)
                        [g1, h1] = onMemory(mu(:,:,:,k), f(:,:,:,k), c, s(:,k), gmu(:,:,:,k,:));
                        g = g + g1;
                        h = h + h1;
                    end
                end
            else
                g = numeric(g);
                if ~par
                    for k=1:nc
                        g1 = onMemory(mu(bbx,bby,bbz,k), f(:,:,:,k), c, s(:,k), gmu(bbx,bby,bbz,k,:));
                        g = g + g1;
                    end
                elseif isa(mu, 'file_array') && isa(gmu, 'file_array') && isa(f, 'file_array')
                    parfor (k=1:nc, par)
                        g1 = onMemory(slicevol(mu, {bbx, bby, bbz, k}), slicevol(f, k, 4), c, s(:,k), slicevol(gmu, {bbx, bby, bbz, k}));
                        g = g + g1;
                    end
                elseif isa(mu, 'file_array') && isa(gmu, 'file_array')
                    parfor (k=1:nc, par)
                        g1 = onMemory(slicevol(mu, {bbx, bby, bbz, k}), f(:,:,:,k), c, s(:,k), slicevol(gmu, {bbx, bby, bbz, k}));
                        g = g + g1;
                    end
                elseif isa(mu, 'file_array') && isa(f, 'file_array')
                    gmu = mu(bbx,bby,bbz,:,:);
                    parfor (k=1:nc, par)
                        g1 = onMemory(slicevol(mu, {bbx, bby, bbz, k}), slicevol(f, k, 4), c, s(:,k), gmu(:,:,:,k,:));
                        g = g + g1;
                    end
                elseif isa(gmu, 'file_array') && isa(f, 'file_array')
                    mu  = mu(bbx,bby,bbz,:);
                    parfor (k=1:nc, par)
                        g1 = onMemory(mu(:,:,:,k), slicevol(f, k, 4), c, s(:,k), slicevol(gmu, {bbx, bby, bbz, k}));
                        g = g + g1;
                    end
                elseif isa(mu, 'file_array')
                    gmu = mu(bbx,bby,bbz,:,:);
                    parfor (k=1:nc, par)
                        g1 = onMemory(slicevol(mu, {bbx, bby, bbz, k}), f(:,:,:,k), c, s(:,k), gmu(:,:,:,k,:));
                        g = g + g1;
                    end
                elseif isa(gmu, 'file_array')
                    mu  = mu(bbx,bby,bbz,:);
                    parfor (k=1:nc, par)
                        g1 = onMemory(mu(:,:,:,k), f(:,:,:,k), c, s(:,k), slicevol(gmu, {bbx, bby, bbz, k}));
                        g = g + g1;
                    end
                elseif isa(f, 'file_array')
                    mu  = mu(bbx,bby,bbz,:);
                    gmu = mu(bbx,bby,bbz,:,:);
                    parfor (k=1:nc, par)
                        g1 = onMemory(mu(:,:,:,k), slicevol(f, k, 4), c, s(:,k), gmu(:,:,:,k,:));
                        g = g + g1;
                    end
                else
                    mu  = mu(bbx,bby,bbz,:);
                    gmu = mu(bbx,bby,bbz,:,:);
                    parfor (k=1:nc, par)
                        g1 = onMemory(mu(:,:,:,k), f(:,:,:,k), c, s(:,k), gmu(:,:,:,k,:));
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
            if hessian
                if ~par
                    for z=1:dlat(3)
                        h1 = onMemory(mu(bbx,bby,oz+z,:), f(:,:,z,:), c(:,:,z), s, [], true);
                        h(:,:,z,:) = h(:,:,z,:) + h1;
                    end
                elseif isa(mu, 'file_array') && isa(f, 'file_array')
                    parfor (z=1:dlat(3), par)
                        h1 = onMemory(slicevol(mu, {bbx,bby,oz+z}), slicevol(f, z, 3), slicevol(c, z, 3), s, [], true);
                        h(:,:,z,:) = h(:,:,z,:) + h1;
                    end
                elseif isa(mu, 'file_array')
                    parfor (z=1:dlat(3), par)
                        h1 = onMemory(slicevol(mu, {bbx,bby,oz+z}), f(:,:,z,:), c(:,:,z), s, [], true);
                        h(:,:,z,:) = h(:,:,z,:) + h1;
                    end
                elseif isa(f, 'file_array')
                    mu = mu(bbx, bby, bbz, :);
                    parfor (z=1:dlat(3), par)
                        h1 = onMemory(mu(:,:,z,:), slicevol(f, z, 3), slicevol(c, z, 3), s, [], true);
                        h(:,:,z,:) = h(:,:,z,:) + h1;
                    end
                else
                    mu = mu(bbx, bby, bbz, :);
                    parfor (z=1:dlat(3), par)
                        h1 = onMemory(mu(:,:,z,:), f(:,:,z,:), c(:,:,z), s, [], true);
                        h(:,:,z,:) = h(:,:,z,:) + h1;
                    end
                end
            elseif nargout > 1
                if ~par
                    for z=1:dlat(3)
                        [g1, h1] = onMemory(mu(bbx,bby,oz+z,:), f(:,:,z,:), c(:,:,z), s);
                        g(:,:,z,:) = g(:,:,z,:) + g1;
                        h(:,:,z,:) = h(:,:,z,:) + h1;
                    end
                elseif isa(mu, 'file_array') && isa(f, 'file_array')
                    parfor (z=1:dlat(3), par)
                        [g1, h1] = onMemory(slicevol(mu, {bbx,bby,oz+z}), slicevol(f, z, 3), slicevol(c, z, 3), s);
                        g(:,:,z,:) = g(:,:,z,:) + g1;
                        h(:,:,z,:) = h(:,:,z,:) + h1;
                    end
                elseif isa(mu, 'file_array')
                    parfor (z=1:dlat(3), par)
                        [g1, h1] = onMemory(slicevol(mu, {bbx,bby,oz+z}), f(:,:,z,:), c(:,:,z), s);
                        g(:,:,z,:) = g(:,:,z,:) + g1;
                        h(:,:,z,:) = h(:,:,z,:) + h1;
                    end
                elseif isa(f, 'file_array')
                    mu = mu(bbx, bby, bbz, :);
                    parfor (z=1:dlat(3), par)
                        [g1, h1] = onMemory(mu(:,:,z,:), slicevol(f, z, 3), slicevol(c, z, 3), s);
                        g(:,:,z,:) = g(:,:,z,:) + g1;
                        h(:,:,z,:) = h(:,:,z,:) + h1;
                    end
                else
                    mu = mu(bbx, bby, bbz, :);
                    parfor (z=1:dlat(3), par)
                        [g1, h1] = onMemory(mu(:,:,z,:), f(:,:,z,:), c(:,:,z), s);
                        g(:,:,z,:) = g(:,:,z,:) + g1;
                        h(:,:,z,:) = h(:,:,z,:) + h1;
                    end
                end
             else
                if ~par
                    for z=1:dlat(3)
                        g1 = onMemory(mu(bbx,bby,oz+z,:), f(:,:,z,:), c(:,:,z), s);
                        g(:,:,z,:) = g(:,:,z,:) + g1;
                    end
                elseif isa(mu, 'file_array') && isa(f, 'file_array')
                    parfor (z=1:dlat(3), par)
                        g1 = onMemory(slicevol(mu, {bbx,bby,oz+z}), slicevol(f, z, 3), slicevol(c, z, 3), s);
                        g(:,:,z,:) = g(:,:,z,:) + g1;
                    end
                elseif isa(mu, 'file_array')
                    parfor (z=1:dlat(3), par)
                        g1 = onMemory(slicevol(mu, {bbx,bby,oz+z}), f(:,:,z,:), c(:,:,z), s);
                        g(:,:,z,:) = g(:,:,z,:) + g1;
                    end
                elseif isa(f, 'file_array')
                    mu = mu(bbx, bby, bbz, :);
                    parfor (z=1:dlat(3), par)
                        g1 = onMemory(mu(:,:,z,:), slicevol(f, z, 3), slicevol(c, z, 3), s);
                        g(:,:,z,:) = g(:,:,z,:) + g1;
                    end
                else
                    mu = mu(bbx, bby, bbz, :);
                    parfor (z=1:dlat(3), par)
                        g1 = onMemory(mu(:,:,z,:), f(:,:,z,:), c(:,:,z), s);
                        g(:,:,z,:) = g(:,:,z,:) + g1;
                    end
                end
            end
        else
            if hessian
                h = numeric(h);
                if ~par
                    for z=1:dlat(3)
                        h1 = onMemory(mu(bbx,bby,oz+z,:), f(:,:,z,:), c(:,:,z), s, gmu(bbx,bby,oz+z,:,:), true);
                        h(:,:,z,:) = h(:,:,z,:) + h1;
                    end
                elseif isa(mu, 'file_array') && isa(gmu, 'file_array') && isa(f, 'file_array')
                    parfor (z=1:dlat(3), par)
                        h1 = onMemory(slicevol(mu, {bbx,bby,oz+z}), slicevol(f, z, 3), slicevol(c, z, 3), s, slicevol(gmu, {bbx,bby,oz+z}), true);
                        h(:,:,z,:) = h(:,:,z,:) + h1;
                    end
                elseif isa(mu, 'file_array') && isa(gmu, 'file_array')
                    parfor (z=1:dlat(3), par)
                        h1 = onMemory(slicevol(mu, {bbx,bby,oz+z}), f(:,:,z,:), c(:,:,z), s, slicevol(gmu, {bbx,bby,oz+z}), true);
                        h(:,:,z,:) = h(:,:,z,:) + h1;
                    end
                elseif isa(mu, 'file_array') && isa(f, 'file_array')
                    gmu = gmu(bbx,bby,bbz,:,:);
                    parfor (z=1:dlat(3), par)
                        h1 = onMemory(slicevol(mu, {bbx,bby,oz+z}), slicevol(f, z, 3), slicevol(c, z, 3), s, gmu(:,:,z,:,:), true);
                        h(:,:,z,:) = h(:,:,z,:) + h1;
                    end
                elseif isa(gmu, 'file_array') && isa(f, 'file_array')
                    mu  = mu(bbx,bby,bbz,:);
                    gmu = gmu(bbx,bby,bbz,:,:);
                    parfor (z=1:dlat(3), par)
                        h1 = onMemory(mu(:,:,z,:), slicevol(f, z, 3), slicevol(c, z, 3), s, slicevol(gmu, {bbx,bby,oz+z}), true);
                        h(:,:,z,:) = h(:,:,z,:) + h1;
                    end
                elseif isa(mu, 'file_array')
                    gmu = gmu(bbx,bby,bbz,:,:);
                    parfor (z=1:dlat(3), par)
                        h1 = onMemory(slicevol(mu, {bbx,bby,oz+z}), f(:,:,z,:), c(:,:,z), s, gmu(:,:,z,:,:), true);
                        h(:,:,z,:) = h(:,:,z,:) + h1;
                    end
                elseif isa(gmu, 'file_array')
                    mu  = mu(bbx,bby,bbz,:);
                    parfor (z=1:dlat(3), par)
                        h1 = onMemory(mu(:,:,z,:), f(:,:,z,:), c(:,:,z), s, slicevol(gmu, {bbx,bby,oz+z}), true);
                        h(:,:,z,:) = h(:,:,z,:) + h1;
                    end
                elseif isa(f, 'file_array')
                    mu  = mu(bbx,bby,bbz,:);
                    gmu = gmu(bbx,bby,bbz,:,:);
                    parfor (z=1:dlat(3), par)
                        h1 = onMemory(mu(:,:,z,:), slicevol(f, z, 3), slicevol(c, z, 3), s, gmu(:,:,z,:,:), true);
                        h(:,:,z,:) = h(:,:,z,:) + h1;
                    end
                else
                    mu  = mu(bbx,bby,bbz,:);
                    gmu = gmu(bbx,bby,bbz,:,:);
                    parfor (z=1:dlat(3), par)
                        h1 = onMemory(mu(:,:,z,:), f(:,:,z,:), c(:,:,z), s, gmu(:,:,z,:,:), true);
                        h(:,:,z,:) = h(:,:,z,:) + h1;
                    end
                end
            elseif nargout > 1
                g = numeric(g);
                h = numeric(h);
                if ~par
                    for z=1:dlat(3)
                        [g1, h1] = onMemory(mu(bbx,bby,oz+z,:), f(:,:,z,:), c(:,:,z), s, gmu(bbx,bby,oz+z,:,:));
                        g(:,:,z,:) = g(:,:,z,:) + g1;
                        h(:,:,z,:) = h(:,:,z,:) + h1;
                    end
                elseif isa(mu, 'file_array') && isa(gmu, 'file_array') && isa(f, 'file_array')
                    parfor (z=1:dlat(3), par)
                        [g1, h1] = onMemory(slicevol(mu, {bbx,bby,oz+z}), slicevol(f, z, 3), slicevol(c, z, 3), s, slicevol(gmu, {bbx,bby,oz+z}));
                        g(:,:,z,:) = g(:,:,z,:) + g1;
                        h(:,:,z,:) = h(:,:,z,:) + h1;
                    end
                elseif isa(mu, 'file_array') && isa(gmu, 'file_array')
                    parfor (z=1:dlat(3), par)
                        [g1, h1] = onMemory(slicevol(mu, {bbx,bby,oz+z}), f(:,:,z,:), c(:,:,z), s, slicevol(gmu, {bbx,bby,oz+z}));
                        g(:,:,z,:) = g(:,:,z,:) + g1;
                        h(:,:,z,:) = h(:,:,z,:) + h1;
                    end
                elseif isa(mu, 'file_array') && isa(f, 'file_array')
                    gmu = gmu(bbx,bby,bbz,:,:);
                    parfor (z=1:dlat(3), par)
                        [g1, h1] = onMemory(slicevol(mu, {bbx,bby,oz+z}), slicevol(f, z, 3), slicevol(c, z, 3), s, gmu(:,:,z,:,:));
                        g(:,:,z,:) = g(:,:,z,:) + g1;
                        h(:,:,z,:) = h(:,:,z,:) + h1;
                    end
                elseif isa(gmu, 'file_array') && isa(f, 'file_array')
                    mu  = mu(bbx,bby,bbz,:);
                    gmu = gmu(bbx,bby,bbz,:,:);
                    parfor (z=1:dlat(3), par)
                        [g1, h1] = onMemory(mu(:,:,z,:), slicevol(f, z, 3), slicevol(c, z, 3), s, slicevol(gmu, {bbx,bby,oz+z}));
                        g(:,:,z,:) = g(:,:,z,:) + g1;
                        h(:,:,z,:) = h(:,:,z,:) + h1;
                    end
                elseif isa(mu, 'file_array')
                    gmu = gmu(bbx,bby,bbz,:,:);
                    parfor (z=1:dlat(3), par)
                        [g1, h1] = onMemory(slicevol(mu, {bbx,bby,oz+z}), f(:,:,z,:), c(:,:,z), s, gmu(:,:,z,:,:));
                        g(:,:,z,:) = g(:,:,z,:) + g1;
                        h(:,:,z,:) = h(:,:,z,:) + h1;
                    end
                elseif isa(gmu, 'file_array')
                    mu  = mu(bbx,bby,bbz,:);
                    parfor (z=1:dlat(3), par)
                        [g1, h1] = onMemory(mu(:,:,z,:), f(:,:,z,:), c(:,:,z), s, slicevol(gmu, {bbx,bby,oz+z}));
                        g(:,:,z,:) = g(:,:,z,:) + g1;
                        h(:,:,z,:) = h(:,:,z,:) + h1;
                    end
                elseif isa(f, 'file_array')
                    mu  = mu(bbx,bby,bbz,:);
                    gmu = gmu(bbx,bby,bbz,:,:);
                    parfor (z=1:dlat(3), par)
                        [g1, h1] = onMemory(mu(:,:,z,:), slicevol(f, z, 3), slicevol(c, z, 3), s, gmu(:,:,z,:,:));
                        g(:,:,z,:) = g(:,:,z,:) + g1;
                        h(:,:,z,:) = h(:,:,z,:) + h1;
                    end
                else
                    mu  = mu(bbx,bby,bbz,:);
                    gmu = gmu(bbx,bby,bbz,:,:);
                    parfor (z=1:dlat(3), par)
                        [g1, h1] = onMemory(mu(:,:,z,:), f(:,:,z,:), c(:,:,z), s, gmu(:,:,z,:,:));
                        g(:,:,z,:) = g(:,:,z,:) + g1;
                        h(:,:,z,:) = h(:,:,z,:) + h1;
                    end
                end
            else
                g = numeric(g);
                if ~par
                    for z=1:dlat(3)
                        g1 = onMemory(mu(bbx,bby,oz+z,:), f(:,:,z,:), c(:,:,z), s, gmu(bbx,bby,oz+z,:,:));
                        g(:,:,z,:) = g(:,:,z,:) + g1;
                    end
                elseif isa(mu, 'file_array') && isa(gmu, 'file_array') && isa(f, 'file_array')
                    parfor (z=1:dlat(3), par)
                        g1 = onMemory(slicevol(mu, {bbx,bby,oz+z}), slicevol(f, z, 3), slicevol(c, z, 3), s, slicevol(gmu, {bbx,bby,oz+z}));
                        g(:,:,z,:) = g(:,:,z,:) + g1;
                    end
                elseif isa(mu, 'file_array') && isa(gmu, 'file_array')
                    parfor (z=1:dlat(3), par)
                        g1 = onMemory(slicevol(mu, {bbx,bby,oz+z}), f(:,:,z,:), c(:,:,z), s, slicevol(gmu, {bbx,bby,oz+z}));
                        g(:,:,z,:) = g(:,:,z,:) + g1;
                    end
                elseif isa(mu, 'file_array') && isa(f, 'file_array')
                    gmu = gmu(bbx,bby,bbz,:,:);
                    parfor (z=1:dlat(3), par)
                        g1 = onMemory(slicevol(mu, {bbx,bby,oz+z}), slicevol(f, z, 3), slicevol(c, z, 3), s, gmu(:,:,z,:,:));
                        g(:,:,z,:) = g(:,:,z,:) + g1;
                    end
                elseif isa(gmu, 'file_array') && isa(f, 'file_array')
                    mu  = mu(bbx,bby,bbz,:);
                    parfor (z=1:dlat(3), par)
                        g1 = onMemory(mu(:,:,z,:), slicevol(f, z, 3), slicevol(c, z, 3), s, slicevol(gmu, {bbx,bby,oz+z}));
                        g(:,:,z,:) = g(:,:,z,:) + g1;
                    end
                elseif isa(mu, 'file_array')
                    gmu = gmu(bbx,bby,bbz,:,:);
                    parfor (z=1:dlat(3), par)
                        g1 = onMemory(slicevol(mu, {bbx,bby,oz+z}), f(:,:,z,:), c(:,:,z), s, gmu(:,:,z,:,:));
                        g(:,:,z,:) = g(:,:,z,:) + g1;
                    end
                elseif isa(gmu, 'file_array')
                    mu  = mu(bbx,bby,bbz,:);
                    parfor (z=1:dlat(3), par)
                        g1 = onMemory(mu(:,:,z,:), f(:,:,z,:), c(:,:,z), s, slicevol(gmu, {bbx,bby,oz+z}));
                        g(:,:,z,:) = g(:,:,z,:) + g1;
                    end
                elseif isa(f, 'file_array')
                    mu  = mu(bbx,bby,bbz,:);
                    gmu = gmu(bbx,bby,bbz,:,:);
                    parfor (z=1:dlat(3), par)
                        g1 = onMemory(mu(:,:,z,:), slicevol(f, z, 3), slicevol(c, z, 3), s, gmu(:,:,z,:,:));
                        g(:,:,z,:) = g(:,:,z,:) + g1;
                    end
                else
                    mu  = mu(bbx,bby,bbz,:);
                    gmu = gmu(bbx,bby,bbz,:,:);
                    parfor (z=1:dlat(3), par)
                        g1 = onMemory(mu(:,:,z,:), f(:,:,z,:), c(:,:,z), s, gmu(:,:,z,:,:));
                        g(:,:,z,:) = g(:,:,z,:) + g1;
                    end
                end
            end
        end
        
    end
    
    % --- Push gradient if needed
    if checkarray(ipsi)
        if hessian || nargin > 1
            h = pushImage(ipsi, h, lat, ...
                'output', h, ...
                'circ',   p.Results.circ, ...
                'loop',   p.Results.loop, ...
                'par',    p.Results.par, ...
                'debug',  p.Results.debug);
            if ~isempty(ga)
                [indv, len] = spm_matcomp('SymIndices', nvec, 'n');
                [hc, h] = deal(h, zeros([lat len], 'single'));
                for z=1:size(h,3)
                    for d=1:nvec
                        for l=d:nvec
                            h(:,:,z,indv(d,l)) = h(:,:,z,indv(d,l)) ...
                                + hc(:,:,z,:) .* ga(:,:,z,:,d) .* ga(:,:,z,:,l);
                        end
                    end
                end
            end
            for z=1:size(h, 3)
                h1 = h(:,:,z,:);
                h1(~isfinite(h1)) = 0;
                h(:,:,z,:) = h1;
            end
        end
        if ~hessian
            g = pushImage(ipsi, g, lat, ...
                'output', g, ...
                'circ',   p.Results.circ, ...
                'loop',   p.Results.loop, ...
                'par',    p.Results.par, ...
                'debug',  p.Results.debug);
            if ~isempty(ga)
                [gc, g] = deal(g, zeros([lat nvec], 'single'));
                for z=1:lat(3)
                    g(:,:,z,:) = -spm_matcomp('Pointwise', ga(:,:,z,:,:), gc(:,:,z,:), 't');
                end
            end
            for z=1:size(g, 3)
                g1 = g(:,:,z,:);
                g1(~isfinite(g1)) = 0;
                g(:,:,z,:) = g1;
            end
        end
    end
    
    % --- Set hessian type
    if isempty(gmu) && isempty(ga)
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
function [g, h] = onMemory(mu, f, c, s, gmu, hessian)
    
    if nargin < 6
        hessian = false;
        if nargin < 5
            gmu = [];
        end
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
    s   = reshape(1./s, [1 1 1 nc]);
    
    if ~hessian
        g  = (bsxfun(@times, c, mu) - f);
        g  = bsxfun(@times, s, g);
        if ~isempty(gmu)
            g = -spm_matcomp('pointwise', gmu, g, 't');
        end
        g(~isfinite(g)) = 0;
    else
        g = [];
    end
    if nargout > 1 || hessian
        h = bsxfun(@times, s, c);
        if numel(c) == 1
            h = repmat(h, lat);
        end
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
    
    if hessian
        [g,h] = deal(h,g);
    end
    
end
