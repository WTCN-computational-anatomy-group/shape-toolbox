function [g, h] = ghMatchingLatent(model, mu, f, gmu, w, varargin)
%__________________________________________________________________________
% 
% Compute gradient/hessian of the **negative** log-likelihood with respect 
% to latent parameters in the principal subspace.
% 
% -------------------------------------------------------------------------
%
% FORMAT [(g), (h)] = ghMatchingLatent(model, mu, f, ga, w, ...)
%
% REQUIRED
% --------
% model - Structure with fields:
%           * 'name'    : 'normal', 'laplace', 'bernoulli' or 'categorical'
%           * ('sigma2'): Normal variance  [1]
%           * ('b')     : Laplace variance [1]
% mu    - Reconstructed probability template [mx my mz nc]
% f     - Observed image [mx my mz nc]
% ga    - Template coefficients spatial gradients
%         > If ipsi provided, [mx my mz nc 3]
%         > Else              [nx ny nz nc 3]
% w     - Principal subspace [nx ny nz 3 nq]
%
% Note: we need either mu warped to image space or f pushed to template
%       space.
%
%  KEYWORD ARGUMENTS
% ------------------
% ipsi    - Inverse (subj to template) warp [mx my mz 3]
%           > If provided, compute "pushed(gradient(pulled))"
% circ    - (Push only) Boundary conditions for pushing [1]
% count   - Pushed voxel count ([nx ny nz nc])
%           > If provided, compute "gradient(pushed)"
% bb      - Bounding box (if different between template and pushed image)
% hessian - Only compute hessian (not gradient) [false]
% loop    - How to split: 'none', 'slice' [auto]
% par     - Parallelise: false/true/number of workers [auto]
% 
% OUTPUT
% ------
% g     - First derivatives w.r.t. latent PC parameters ([nq])
% h     - Second derivatives w.r.t. latent PC parameters ([nq nq])
%__________________________________________________________________________

    % --- Parse inputs
    p = inputParser;
    p.FunctionName = 'ghMatchingLatent';
    p.addRequired('model',  @(X) isstruct(X) && isfield(X, 'name'));
    p.addRequired('mu',     @checkarray);
    p.addRequired('f',      @checkarray);
    p.addRequired('gmu',    @checkarray);
    p.addRequired('w',      @checkarray);
    p.addParameter('ipsi',     [],      @checkarray);
    p.addParameter('circ',    1);
    p.addParameter('count',    [],      @checkarray);
    p.addParameter('bb',       struct,  @isstruct);
    p.addParameter('hessian',  false,   @islogical);
    p.addParameter('loop',     '',      @ischar);
    p.addParameter('par',      true,    @isscalar);
    p.addParameter('debug',    false,   @isscalar);
    p.parse(model, mu, f, gmu, w, varargin{:});
    loop     = p.Results.loop;
    par      = p.Results.par;
    debug    = p.Results.debug;
    bb       = p.Results.bb;
    c        = p.Results.count;
    ipsi     = p.Results.ipsi;
    circ     = p.Results.circ;
    hessian  = p.Results.hessian;
    
    if p.Results.debug, fprintf('* ghMatchingLatent\n'); end
    
    % --- Optimise parallelisation and splitting schemes
    [par, loop] = autoParLoop(par, loop, isa(w, 'file_array'), size(w,3));
    
    % --- Default bounding box
    dim = [size(mu) 1];
    if ~isfield(bb, 'x')
        bb.x = 1:dim(1);
    end
    if ~isfield(bb, 'y')
        bb.y = 1:dim(2);
    end
    if ~isfield(bb, 'z')
        bb.z = 1:dim(3);
    end
    bbx = bb.x;
    bby = bb.y;
    bbz = bb.z;
    oz = bbz(1) - 1;

    % --- Allocate arrays for grad and hessian
    do_gradient = ~hessian;
    do_hessian = nargout > 1 || hessian;
    nk  = size(w,5);
    if do_gradient
        g = zeros([nk 1], 'double');
    end
    if do_hessian
        h = zeros(nk, 'double');
    end
    if do_hessian && do_gradient
        if isempty(ipsi)
            nout = 2;
        else
            nout = 3;
        end
    else
        if isempty(ipsi)
            nout = 1;
        else
            nout = 2;
        end
    end

    % ---------------------------------------------------------------------
    %   Case pulled gradient : Pushed(Gradient(Pulled)) 
    % ---------------------------------------------------------------------
    if ~isempty(ipsi)
        
        % Compute gradient/hessian w.r.t. template coefficients
        % -----------------------------------------------------
        
        output = {};
        if isa(mu, 'file_array')
            [path, name, ext] = fileparts(mu.fname);
            if do_gradient
                copyfile(mu.fname, fullfile(path, [name '_gca' ext]));
                gca = nifti(fullfile(path, [name '_gca' ext]));
                gca = gca.dat;
                output = [output {gca}];
            end
            if do_hessian
                copyfile(mu.fname, fullfile(path, [name '_hca' ext]));
                hca = nifti(fullfile(path, [name '_hca' ext]));
                hca = hca.dat;
                output = [output {hca}];
            end
        end
        wlat = [size(w) 1];
        wlat = wlat(1:3);
        [argout{1:nout}] = ghMatchingVel(model, mu, f, ...
            'ipsi',    ipsi, ...
            'lat',     wlat, ...
            'circ',    circ, ...
            'bb',      bb, ...
            'hessian', ~do_gradient, ...
            'loop',    loop, ...
            'par',     par, ...
            'output',  output);
        if do_gradient
            gca = argout{1};
            if do_hessian
                hca = argout{2};
                htype = argout{3};
                nout = 2;
            else
                hca = zeros([0 0 wlat(3) 0]); % sliceable empty array
                htype = '';
                nout = 1;
            end
        elseif do_hessian
            hca = argout{1};
            htype = argout{2};
            gca = zeros([0 0 wlat(3) 0]); % sliceable empty array
            nout = 1;
        end
        clear model mu f 
        
        % Multiplication with PC*ga
        % -------------------------
        
        % --- If loop on slices
        if strcmpi(loop, 'slice')
            nz = wlat(3);
            if debug
                if par
                    fprintf('  - Parallelise on slices\n')
                else
                    fprintf('  - Serialise on slices\n')
                end
            end
            fa_gca = isa(gca, 'file_array');
            fa_hca = isa(hca, 'file_array');
            fa_gmu = isa(gmu, 'file_array');
            fa_w   = isa(w,   'file_array');
            
            if ~par
                for z=1:nz
                    [argout{1:nout}] = onMemoryPulled( ...
                        gca(:,:,z,:), ...
                        hca(:,:,z,:), htype, ...
                        gmu(:,:,z,:,:), ...
                        w(:,:,z,:,:), ...
                        'loop',     'none', ...
                        'hessian',  ~do_gradient, ...
                        'par',      false);
                    if do_gradient
                        g = g + argout{1};
                        if do_hessian
                            h = h + argout{2};
                        end
                    elseif do_hessian
                            h = h + argout{1};
                    end
                end
            elseif  fa_gca &&  fa_hca &&  fa_gmu &&  fa_w
                parfor (z=1:nz, par)
                    argout = {};
                    [argout{1:nout}] = onMemoryPulled( ...
                        slicevol(gca, z, 3), ...
                        slicevol(hca, z, 3), htype, ...
                        slicevol(gmu, z, 3), ...
                        slicevol(w, z, 3), ...
                        'loop',     'none', ...
                        'hessian',  ~do_gradient, ...
                        'par',      false);
                    if do_gradient
                        g = g + argout{1};
                        if do_hessian
                            h = h + argout{2};
                        end
                    elseif do_hessian
                            h = h + argout{1};
                    end
                end
            elseif  fa_gca &&  fa_hca &&  fa_gmu && ~fa_w
                parfor (z=1:nz, par)
                    argout = {};
                    [argout{1:nout}] = onMemoryPulled( ...
                        slicevol(gca, z, 3), ...
                        slicevol(hca, z, 3), htype, ...
                        slicevol(gmu, z, 3), ...
                        w(:,:,z,:,:), ...
                        'loop',     'none', ...
                        'hessian',  ~do_gradient, ...
                        'par',      false);
                    if do_gradient
                        g = g + argout{1};
                        if do_hessian
                            h = h + argout{2};
                        end
                    elseif do_hessian
                            h = h + argout{1};
                    end
                end
            elseif  fa_gca &&  fa_hca && ~fa_gmu &&  fa_w
                parfor (z=1:nz, par)
                    argout = {};
                    [argout{1:nout}] = onMemoryPulled( ...
                        slicevol(gca, z, 3), ...
                        slicevol(hca, z, 3), htype, ...
                        gmu(:,:,z,:,:), ...
                        slievol(w, z, 3), ...
                        'loop',     'none', ...
                        'hessian',  ~do_gradient, ...
                        'par',      false);
                    if do_gradient
                        g = g + argout{1};
                        if do_hessian
                            h = h + argout{2};
                        end
                    elseif do_hessian
                            h = h + argout{1};
                    end
                end
            elseif  fa_gca && ~fa_hca &&  fa_gmu &&  fa_w
                parfor (z=1:nz, par)
                    argout = {};
                    [argout{1:nout}] = onMemoryPulled( ...
                        slicevol(gca, z, 3), ...
                        hca(:,:,z,:), htype, ...
                        slicevol(gmu, z, 3), ...
                        slicevol(w, z, 3), ...
                        'loop',     'none', ...
                        'hessian',  ~do_gradient, ...
                        'par',      false);
                    if do_gradient
                        g = g + argout{1};
                        if do_hessian
                            h = h + argout{2};
                        end
                    elseif do_hessian
                            h = h + argout{1};
                    end
                end
            elseif ~fa_gca &&  fa_hca &&  fa_gmu &&  fa_w
                parfor (z=1:nz, par)
                    argout = {};
                    [argout{1:nout}] = onMemoryPulled( ...
                        gca(:,:,z,:), ...
                        slicevol(hca, 3), htype, ...
                        slicevol(gmu, z, 3), ...
                        slicevol(w, z, 3), ...
                        'loop',     'none', ...
                        'hessian',  ~do_gradient, ...
                        'par',      false);
                    if do_gradient
                        g = g + argout{1};
                        if do_hessian
                            h = h + argout{2};
                        end
                    elseif do_hessian
                            h = h + argout{1};
                    end
                end
            elseif  fa_gca &&  fa_hca && ~fa_gmu && ~fa_w
                parfor (z=1:nz, par)
                    argout = {};
                    [argout{1:nout}] = onMemoryPulled( ...
                        slicevol(gca, z, 3), ...
                        slicevol(hca, z, 3), htype, ...
                        gmu(:,:,z,:,:), ...
                        w(:,:,z,:,:), ...
                        'loop',     'none', ...
                        'hessian',  ~do_gradient, ...
                        'par',      false);
                    if do_gradient
                        g = g + argout{1};
                        if do_hessian
                            h = h + argout{2};
                        end
                    elseif do_hessian
                            h = h + argout{1};
                    end
                end
            elseif  fa_gca && ~fa_hca &&  fa_gmu && ~fa_w
                parfor (z=1:nz, par)
                    argout = {};
                    [argout{1:nout}] = onMemoryPulled( ...
                        slicevol(gca, z, 3), ...
                        hca(:,:,z,:), htype, ...
                        slicevol(gmu, z, 3), ...
                        w(:,:,z,:,:), ...
                        'loop',     'none', ...
                        'hessian',  ~do_gradient, ...
                        'par',      false);
                    if do_gradient
                        g = g + argout{1};
                        if do_hessian
                            h = h + argout{2};
                        end
                    elseif do_hessian
                            h = h + argout{1};
                    end
                end
            elseif ~fa_gca &&  fa_hca &&  fa_gmu && ~fa_w
                parfor (z=1:nz, par)
                    argout = {};
                    [argout{1:nout}] = onMemoryPulled( ...
                        gca(:,:,z,:), ...
                        slice(hca, z, 3), htype, ...
                        slicevol(gmu, z, 3), ...
                        w(:,:,z,:,:), ...
                        'loop',     'none', ...
                        'hessian',  ~do_gradient, ...
                        'par',      false);
                    if do_gradient
                        g = g + argout{1};
                        if do_hessian
                            h = h + argout{2};
                        end
                    elseif do_hessian
                            h = h + argout{1};
                    end
                end
            elseif  fa_gca && ~fa_hca && ~fa_gmu &&  fa_w
                parfor (z=1:nz, par)
                    argout = {};
                    [argout{1:nout}] = onMemoryPulled( ...
                        slicevol(gca, z, 3), ...
                        hca(:,:,z,:), htype, ...
                        gmu(:,:,z,:,:), ...
                        slicevol(w, z, 3), ...
                        'loop',     'none', ...
                        'hessian',  ~do_gradient, ...
                        'par',      false);
                    if do_gradient
                        g = g + argout{1};
                        if do_hessian
                            h = h + argout{2};
                        end
                    elseif do_hessian
                            h = h + argout{1};
                    end
                end
            elseif ~fa_gca &&  fa_hca && ~fa_gmu &&  fa_w
                parfor (z=1:nz, par)
                    argout = {};
                    [argout{1:nout}] = onMemoryPulled( ...
                        gca(:,:,z,:), ...
                        slicevol(hca, z, 3), htype, ...
                        gmu(:,:,z,:,:), ...
                        slicevol(w, z, 3), ...
                        'loop',     'none', ...
                        'hessian',  ~do_gradient, ...
                        'par',      false);
                    if do_gradient
                        g = g + argout{1};
                        if do_hessian
                            h = h + argout{2};
                        end
                    elseif do_hessian
                            h = h + argout{1};
                    end
                end
            elseif ~fa_gca && ~fa_hca &&  fa_gmu &&  fa_w
                parfor (z=1:nz, par)
                    argout = {};
                    [argout{1:nout}] = onMemoryPulled( ...
                        gca(:,:,z,:), ...
                        hca(:,:,z,:), htype, ...
                        slicevol(gmu, z, 3), ...
                        slicevol(w, z, 3), ...
                        'loop',     'none', ...
                        'hessian',  ~do_gradient, ...
                        'par',      false);
                    if do_gradient
                        g = g + argout{1};
                        if do_hessian
                            h = h + argout{2};
                        end
                    elseif do_hessian
                            h = h + argout{1};
                    end
                end
            elseif  fa_gca && ~fa_hca && ~fa_gmu && ~fa_w
                parfor (z=1:nz, par)
                    argout = {};
                    [argout{1:nout}] = onMemoryPulled( ...
                        slicevol(gca, z, 3), ...
                        hca(:,:,z,:), htype, ...
                        gmu(:,:,z,:,:), ...
                        slicevol(w, z, 3), ...
                        'loop',     'none', ...
                        'hessian',  ~do_gradient, ...
                        'par',      false);
                    if do_gradient
                        g = g + argout{1};
                        if do_hessian
                            h = h + argout{2};
                        end
                    elseif do_hessian
                            h = h + argout{1};
                    end
                end
            elseif ~fa_gca &&  fa_hca && ~fa_gmu && ~fa_w
                parfor (z=1:nz, par)
                    argout = {};
                    [argout{1:nout}] = onMemoryPulled( ...
                        gca(:,:,z,:), ...
                        slicevol(hca, z, 3), htype, ...
                        gmu(:,:,z,:,:), ...
                        w(:,:,z,:,:), ...
                        'loop',     'none', ...
                        'hessian',  ~do_gradient, ...
                        'par',      false);
                    if do_gradient
                        g = g + argout{1};
                        if do_hessian
                            h = h + argout{2};
                        end
                    elseif do_hessian
                            h = h + argout{1};
                    end
                end
            elseif ~fa_gca && ~fa_hca &&  fa_gmu && ~fa_w
                parfor (z=1:nz, par)
                    argout = {};
                    [argout{1:nout}] = onMemoryPulled( ...
                        gca(:,:,z,:), ...
                        hca(:,:,z,:), htype, ...
                        slicevol(gmu, z, 3), ...
                        w(:,:,z,:,:), ...
                        'loop',     'none', ...
                        'hessian',  ~do_gradient, ...
                        'par',      false);
                    if do_gradient
                        g = g + argout{1};
                        if do_hessian
                            h = h + argout{2};
                        end
                    elseif do_hessian
                            h = h + argout{1};
                    end
                end
            elseif ~fa_gca && ~fa_hca && ~fa_gmu &&  fa_w
                parfor (z=1:nz, par)
                    argout = {};
                    [argout{1:nout}] = onMemoryPulled( ...
                        gca(:,:,z,:), ...
                        hca(:,:,z,:), htype, ...
                        gmu(:,:,z,:,:), ...
                        slicevol(w, z, 3), ...
                        'loop',     'none', ...
                        'hessian',  ~do_gradient, ...
                        'par',      false);
                    if do_gradient
                        g = g + argout{1};
                        if do_hessian
                            h = h + argout{2};
                        end
                    elseif do_hessian
                            h = h + argout{1};
                    end
                end
            else % ~fa_gca && ~fa_hca && ~fa_gmu && ~fa_w
                parfor (z=1:nz, par)
                    argout = {};
                    [argout{1:nout}] = onMemoryPulled( ...
                        gca(:,:,z,:), ...
                        hca(:,:,z,:), htype, ...
                        gmu(:,:,z,:,:), ...
                        w(:,:,z,:,:), ...
                        'loop',     'none', ...
                        'hessian',  ~do_gradient, ...
                        'par',      false);
                    if do_gradient
                        g = g + argout{1};
                        if do_hessian
                            h = h + argout{2};
                        end
                    elseif do_hessian
                            h = h + argout{1};
                    end
                end
            end
        % --- No loop
        else
            if debug
                fprintf('  - No loop\n')
            end
            [argout{1:nout}] = onMemoryPulled( ...
                gca, hca, htype, gmu, w, ...
                'loop',     loop, ...
                'hessian',  ~do_gradient, ...
                'par',      par, ...
                'debug',    debug);
            if do_gradient
                g = argout{1};
                if do_hessian
                    h = argout{2};
                end
            elseif do_hessian
                    h = argout{1};
            end
        end
          
        if isa(gca, 'file_array')
            rmarray(gca);
        end
        if isa(hca, 'file_array')
            rmarray(hca);
        end
        clear gca hca
    % ---------------------------------------------------------------------
    %   Case pushed gradient : Gradient(Pushed)
    % ---------------------------------------------------------------------
    else

        % --- If loop on slices
        if strcmpi(loop, 'slice')
            nz = length(bbz);
            if debug
                if par
                    fprintf('  - Parallelise on slices\n')
                else
                    fprintf('  - Serialise on slices\n')
                end
            end
            % Parfor has stupid rules:
            % 1) it loads all slices in memory at once, so file_array must be
            %    sliced with a subfunction to avoid memory explosion.
            % 2) it only sends subarrays to the workers if they specificaly
            %    appear in the form a(:,:,z,:) directly in the loop body 
            %    (we cannot use stored offsets).
            % This forces me to write a different loop for all different cases
            % of "on memory" vs "file array".
            fa_mu  = isa(mu,  'file_array');
            fa_f   = isa(f,   'file_array');
            fa_gmu = isa(gmu, 'file_array');
            fa_w   = isa(w,   'file_array');
            % ------------------------
            %   GRADIENT AND HESSIAN
            % ------------------------
            if ~par
                for z=1:nz
                    [argout{1:nout}] = onMemoryPushed(model, ...
                        mu(bbx,bby,oz+z,:), ...
                        f(:,:,z,:), c(:,:,z), ...
                        gmu(bbx,bby,oz+z,:,:), ...
                        w(bbx,bby,oz+z,:,:), ...
                        'loop',     'none', ...
                        'hessian',  ~do_gradient, ...
                        'par',      false);
                    if do_gradient
                        g = g + argout{1};
                        if do_hessian
                            h = h + argout{2};
                        end
                    elseif do_hessian
                            h = h + argout{1};
                    end
                end
            elseif  fa_mu &&  fa_f &&  fa_gmu &&  fa_w
                parfor (z=1:nz, par)
                    argout = {};
                    [argout{1:nout}] = onMemoryPushed(model, ...
                        slicevol(mu, {bbx,bby,oz+z}), ...
                        slicevol(f, z, 3), slicevol(c, z, 3), ...
                        slicevol(gmu, {bbx,bby,oz+z}), ...
                        slicevol(w, {bbx,bby,oz+z}), ...
                        'loop',     'none', ...
                        'hessian',  ~do_gradient, ...
                        'par',      false);
                    if do_gradient
                        g = g + argout{1};
                        if do_hessian
                            h = h + argout{2};
                        end
                    elseif do_hessian
                            h = h + argout{1};
                    end
                end
            elseif  fa_mu &&  fa_f &&  fa_gmu && ~fa_w
                w   = w(bbx,bby,bbz,:,:);
                parfor (z=1:nz, par)
                    argout = {};
                    [argout{1:nout}] = onMemoryPushed(model, ...
                        slicevol(mu, {bbx, bby, oz+z}), ...
                        slicevol(f, z, 3), slicevol(c, z, 3), ...
                        slicevol(gmu, {bbx, bby, oz+z}), ...
                        w(:,:,z,:,:), ...
                        'loop',     'none', ...
                        'hessian',  ~do_gradient, ...
                        'par',      false);
                    if do_gradient
                        g = g + argout{1};
                        if do_hessian
                            h = h + argout{2};
                        end
                    elseif do_hessian
                            h = h + argout{1};
                    end
                end
            elseif  fa_mu &&  fa_f && ~fa_gmu &&  fa_w
                gmu = gmu(bbx,bby,bbz,:,:);
                parfor (z=1:nz, par)
                    argout = {};
                    [argout{1:nout}] = onMemoryPushed(model, ...
                        slicevol(mu, {bbx, bby, oz+z}), ...
                        slicevol(f, z, 3), slicevol(c, z, 3), ...
                        gmu(:,:,z,:,:), ...
                        slievol(w, {bbx, bby, oz+z}), ...
                        'loop',     'none', ...
                        'hessian',  ~do_gradient, ...
                        'par',      false);
                    if do_gradient
                        g = g + argout{1};
                        if do_hessian
                            h = h + argout{2};
                        end
                    elseif do_hessian
                            h = h + argout{1};
                    end
                end
            elseif  fa_mu && ~fa_f &&  fa_gmu &&  fa_w
                parfor (z=1:nz, par)
                    argout = {};
                    [argout{1:nout}] = onMemoryPushed(model, ...
                        slicevol(mu, {bbx, bby, oz+z}), ...
                        f(:,:,z,:), c(:,:,z), ...
                        slicevol(gmu, {bbx, bby, oz+z}), ...
                        slicevol(w, {bbx, bby, oz+z}), ...
                        'loop',     'none', ...
                        'hessian',  ~do_gradient, ...
                        'par',      false);
                    if do_gradient
                        g = g + argout{1};
                        if do_hessian
                            h = h + argout{2};
                        end
                    elseif do_hessian
                            h = h + argout{1};
                    end
                end
            elseif ~fa_mu &&  fa_f &&  fa_gmu &&  fa_w
                mu  = mu(bbx,bby,bbz,:);
                parfor (z=1:nz, par)
                    argout = {};
                    [argout{1:nout}] = onMemoryPushed(model, ...
                        mu(:,:,z,:), ...
                        slicevol(fz, 3), slicevol(c, z, 3), ...
                        slicevol(gmu, {bbx, bby, oz+z}), ...
                        slicevol(w, {bbx, bby, oz+z}), ...
                        'loop',     'none', ...
                        'hessian',  ~do_gradient, ...
                        'par',      false);
                    if do_gradient
                        g = g + argout{1};
                        if do_hessian
                            h = h + argout{2};
                        end
                    elseif do_hessian
                            h = h + argout{1};
                    end
                end
            elseif  fa_mu &&  fa_f && ~fa_gmu && ~fa_w
                gmu = gmu(bbx,bby,bbz,:,:);
                w   = w(bbx,bby,bbz,:,:);
                parfor (z=1:nz, par)
                    argout = {};
                    [argout{1:nout}] = onMemoryPushed(model, ...
                        slicevol(mu, {bbx, bby, oz+z}), ...
                        slicevol(f, z, 3), slicevol(c, z, 3), ...
                        gmu(:,:,z,:,:), ...
                        w(:,:,z,:,:), ...
                        'loop',     'none', ...
                        'hessian',  ~do_gradient, ...
                        'par',      false);
                    if do_gradient
                        g = g + argout{1};
                        if do_hessian
                            h = h + argout{2};
                        end
                    elseif do_hessian
                            h = h + argout{1};
                    end
                end
            elseif  fa_mu && ~fa_f &&  fa_gmu && ~fa_w
                w   = w(bbx,bby,bbz,:,:);
                parfor (z=1:nz, par)
                    argout = {};
                    [argout{1:nout}] = onMemoryPushed(model, ...
                        slicevol(mu, {bbx, bby, oz+z}), ...
                        f(:,:,z,:), c(:,:,z), ...
                        slicevol(gmu, {bbx, bby, oz+z}), ...
                        w(:,:,z,:,:), ...
                        'loop',     'none', ...
                        'hessian',  ~do_gradient, ...
                        'par',      false);
                    if do_gradient
                        g = g + argout{1};
                        if do_hessian
                            h = h + argout{2};
                        end
                    elseif do_hessian
                            h = h + argout{1};
                    end
                end
            elseif ~fa_mu &&  fa_f &&  fa_gmu && ~fa_w
                mu  = mu(bbx,bby,bbz,:);
                w   = w(bbx,bby,bbz,:,:);
                parfor (z=1:nz, par)
                    argout = {};
                    [argout{1:nout}] = onMemoryPushed(model, ...
                        mu(:,:,z,:), ...
                        slice(f, z, 3), slicevol(c, z, 3), ...
                        slicevol(gmu, {bbx, bby, oz+z}), ...
                        w(:,:,z,:,:), ...
                        'loop',     'none', ...
                        'hessian',  ~do_gradient, ...
                        'par',      false);
                    if do_gradient
                        g = g + argout{1};
                        if do_hessian
                            h = h + argout{2};
                        end
                    elseif do_hessian
                            h = h + argout{1};
                    end
                end
            elseif  fa_mu && ~fa_f && ~fa_gmu &&  fa_w
                gmu = gmu(bbx,bby,bbz,:,:);
                parfor (z=1:nz, par)
                    argout = {};
                    [argout{1:nout}] = onMemoryPushed(model, ...
                        slicevol(mu, {bbx, bby, oz+z}), ...
                        f(:,:,z,:), c(:,:,z), ...
                        gmu(:,:,z,:,:), ...
                        slicevol(w, {bbx, bby, oz+z}), ...
                        'loop',     'none', ...
                        'hessian',  ~do_gradient, ...
                        'par',      false);
                    if do_gradient
                        g = g + argout{1};
                        if do_hessian
                            h = h + argout{2};
                        end
                    elseif do_hessian
                            h = h + argout{1};
                    end
                end
            elseif ~fa_mu &&  fa_f && ~fa_gmu &&  fa_w
                mu  = mu(bbx,bby,bbz,:);
                gmu = gmu(bbx,bby,bbz,:,:);
                parfor (z=1:nz, par)
                    argout = {};
                    [argout{1:nout}] = onMemoryPushed(model, ...
                        mu(:,:,z,:), ...
                        slicevol(f, z, 3), slicevol(c, z, 3), ...
                        gmu(:,:,z,:,:), ...
                        slicevol(w, {bbx, bby, oz+z}), ...
                        'loop',     'none', ...
                        'hessian',  ~do_gradient, ...
                        'par',      false);
                    if do_gradient
                        g = g + argout{1};
                        if do_hessian
                            h = h + argout{2};
                        end
                    elseif do_hessian
                            h = h + argout{1};
                    end
                end
            elseif ~fa_mu && ~fa_f &&  fa_gmu &&  fa_w
                mu  = mu(bbx,bby,bbz,:);
                parfor (z=1:nz, par)
                    argout = {};
                    [argout{1:nout}] = onMemoryPushed(model, ...
                        mu(:,:,z,:), ...
                        f(:,:,z,:), c(:,:,z), ...
                        slicevol(gmu, {bbx, bby, oz+z}), ...
                        slicevol(w, {bbx, bby, oz+z}), ...
                        'loop',     'none', ...
                        'hessian',  ~do_gradient, ...
                        'par',      false);
                    if do_gradient
                        g = g + argout{1};
                        if do_hessian
                            h = h + argout{2};
                        end
                    elseif do_hessian
                            h = h + argout{1};
                    end
                end
            elseif  fa_mu && ~fa_f && ~fa_gmu && ~fa_w
                gmu = gmu(bbx,bby,bbz,:,:);
                parfor (z=1:nz, par)
                    argout = {};
                    [argout{1:nout}] = onMemoryPushed(model, ...
                        slicevol(mu, {bbx, bby, oz+z}), ...
                        f(:,:,z,:), c(:,:,z), ...
                        gmu(:,:,z,:,:), ...
                        slicevol(w, {bbx, bby, oz+z}), ...
                        'loop',     'none', ...
                        'hessian',  ~do_gradient, ...
                        'par',      false);
                    if do_gradient
                        g = g + argout{1};
                        if do_hessian
                            h = h + argout{2};
                        end
                    elseif do_hessian
                            h = h + argout{1};
                    end
                end
            elseif ~fa_mu &&  fa_f && ~fa_gmu && ~fa_w
                mu  = mu(bbx,bby,bbz,:);
                gmu = gmu(bbx,bby,bbz,:,:);
                w   = w(bbx,bby,bbz,:,:);
                parfor (z=1:nz, par)
                    argout = {};
                    [argout{1:nout}] = onMemoryPushed(model, ...
                        mu(:,:,z,:), ...
                        slicevol(f, z, 3), slicevol(c, z, 3), ...
                        gmu(:,:,z,:,:), ...
                        w(:,:,z,:,:), ...
                        'loop',     'none', ...
                        'hessian',  ~do_gradient, ...
                        'par',      false);
                    if do_gradient
                        g = g + argout{1};
                        if do_hessian
                            h = h + argout{2};
                        end
                    elseif do_hessian
                            h = h + argout{1};
                    end
                end
            elseif ~fa_mu && ~fa_f &&  fa_gmu && ~fa_w
                mu  = mu(bbx,bby,bbz,:);
                w   = w(bbx,bby,bbz,:,:);
                parfor (z=1:nz, par)
                    argout = {};
                    [argout{1:nout}] = onMemoryPushed(model, ...
                        mu(:,:,z,:), ...
                        f(:,:,z,:), c(:,:,z), ...
                        slicevol(gmu, {bbx, bby, oz+z}), ...
                        w(:,:,z,:,:), ...
                        'loop',     'none', ...
                        'hessian',  ~do_gradient, ...
                        'par',      false);
                    if do_gradient
                        g = g + argout{1};
                        if do_hessian
                            h = h + argout{2};
                        end
                    elseif do_hessian
                            h = h + argout{1};
                    end
                end
            elseif ~fa_mu && ~fa_f && ~fa_gmu &&  fa_w
                mu  = mu(bbx,bby,bbz,:);
                gmu = gmu(bbx,bby,bbz,:,:);
                parfor (z=1:nz, par)
                    argout = {};
                    [argout{1:nout}] = onMemoryPushed(model, ...
                        mu(:,:,z,:), ...
                        f(:,:,z,:), c(:,:,z), ...
                        gmu(:,:,z,:,:), ...
                        slicevol(w, {bbx, bby, oz+z}), ...
                        'loop',     'none', ...
                        'hessian',  ~do_gradient, ...
                        'par',      false);
                    if do_gradient
                        g = g + argout{1};
                        if do_hessian
                            h = h + argout{2};
                        end
                    elseif do_hessian
                            h = h + argout{1};
                    end
                end
            else % ~fa_mu && ~fa_f && ~fa_gmu && ~fa_w
                mu  = mu(bbx,bby,bbz,:);
                gmu = gmu(bbx,bby,bbz,:,:);
                w   = w(bbx,bby,bbz,:,:);
                parfor (z=1:nz, par)
                    argout = {};
                    [argout{1:nout}] = onMemoryPushed(model, ...
                        mu(:,:,z,:), ...
                        f(:,:,z,:), c(:,:,z), ...
                        gmu(:,:,z,:,:), ...
                        w(:,:,z,:,:), ...
                        'loop',     'none', ...
                        'hessian',  ~do_gradient, ...
                        'par',      false);
                    if do_gradient
                        g = g + argout{1};
                        if do_hessian
                            h = h + argout{2};
                        end
                    elseif do_hessian
                            h = h + argout{1};
                    end
                end
            end
        % --- No loop
        else
            if debug
                fprintf('  - No loop\n')
            end
            mu  = mu(bbx,bby,bbz,:);
            gmu = gmu(bbx,bby,bbz,:,:);
            w   = w(bbx,bby,bbz,:,:);
            [argout{1:nargout}] = onMemoryPushed(model, mu, f, c, gmu, w, ...
                'loop',     loop, ...
                'hessian',  ~do_gradient, ...
                'par',      par, ...
                'debug',    debug);
            if do_gradient
                g = argout{1};
                if do_hessian
                    h = argout{2};
                end
            elseif do_hessian
                    h = argout{1};
            end
        end
        
    end % < switch cases
    
    clear argout

    if do_gradient
        g(~isfinite(g)) = 0;
    end
    if do_hessian
        % Insure symmetric hessian
        h = (h + h')/2;
        h(~isfinite(h)) = 0;
    end
    
    if hessian
        g = [];
        [g, h] = deal(h, g);
    end
end

function [g, h] = onMemoryPushed(model, mu, f, c, gmu, w, varargin)


    % --- Multiply PCs with template gradients
    w = -spm_matcomp('Pointwise', single(numeric(gmu)), single(numeric(w)));
    clear gmu
    % => size(w) = [nx ny nz nclasses nlatent]

    % --- Compute grad/hess w.r.t. the complete velocity
    % (fast version that does not perform grad multiplication)
    hessian = find(strcmpi('hessian', varargin));
    hessian = ~isempty(hessian) && varargin{hessian+1};
    if nargout == 1
        if hessian
            [h, htype] = ghMatchingVel(model, mu, f, 'count', c, varargin{:});
        else
            g = ghMatchingVel(model, mu, f, 'count', c, varargin{:});
        end
    else
        [g, h, htype] = ghMatchingVel(model, mu, f, 'count', c, varargin{:});
    end
    clear mu f c

    % --- Compute grad/hess w.r.t. the latent (low-dim) coordinates
    if nargout == 1
        if hessian
            h = vel2latGradHessMatching(w, [], h, htype);
        else
            g = vel2latGradHessMatching(w, g);
        end
    else
        [g, h] = vel2latGradHessMatching(w, g, h, htype);
    end
    clear w

    if hessian
        g = [];
        [g, h] = deal(h, g);
    end
end


function [g, h] = onMemoryPulled(g, h, htype, ga, w, varargin)
% Compute slice-wise gradient and hessian when the deriatives w.r.t.
% template parameters were computed in subject space (pulled) and then 
% pushed to template space.
%
% g     - Pushed gradient w.r.t. changes in warped template [nx ny nz k]
% h     - Pushed Hessian w.r.t. changes in warped template [nx ny nz k k]
% htype - 'symtensor' or 'diagonal'
% ga    - Template coefficients spatial gradients, in template space
%         [nx ny nz k 3]
% w     - Principal subspace [nx ny nz 3 nq]

    % --- Multiply PCs with template gradients
    w = -spm_matcomp('Pointwise', single(numeric(ga)), single(numeric(w)));
    clear ga
    % => size(w) = [nx ny nz nclasses nlatent]


    % --- Compute grad/hess w.r.t. the latent (low-dim) coordinates
    hessian = find(strcmpi('hessian', varargin));
    hessian = ~isempty(hessian) && varargin{hessian+1};
    if nargout == 1
        if hessian
            h = vel2latGradHessMatching(w, [], h, htype);
        else
            g = vel2latGradHessMatching(w, g);
        end
    else
        [g, h] = vel2latGradHessMatching(w, g, h, htype);
    end
    clear w

    if hessian
        g = [];
        [g, h] = deal(h, g);
    end
end

function [g, h] = vel2latGradHessMatching(w, g, h1, htype)
% FORMAT [(g), (h)] = vel2latGradHessMatching(w, g, h, htype)
% w     - Shape subspace basis (pre-multiplied by the template gradients)
% g     - Gradient w.r.t. the full velocity
% h     - Hessian w.r.t. the full velocity
% htype - Type of the hessian approximation w.r.t. classes:
%         'diagonal' or 'symtensor'
%
% Apply the chain rule to convert Grad/Hess w.r.t. the full velocity to
% Grad/Hes w.r.t. the latent coordinates.

    % --- Dim info
    dim         = [size(w) 1 1 1];
    dim         = dim(1:5);
    dim_lattice = dim(1:3);
    dim_classes = dim(4);   % Number of classes (Pre-multiplied W(xi,k) * -Gmu(xi))
    dim_latent  = dim(5);
    
    % --- Default arguments
    do_hessian  = nargin > 2;
    do_gradient = ~isempty(g);
    if do_hessian && nargin < 4
        % Try to guess htype
        if issame(size(h1, 4), dim_latent)
            htype = 'diagonal';
        else
            htype = 'symtensor';
        end
    end
    
    % --- Gradient
    % G = w' * g1
    w = reshape(w, [], dim_latent); % Matrix form
    if do_gradient
        g = double(w' * g(:));
    end
    
    % --- Hessian
    if do_hessian
        h = zeros(dim_latent, 'double');
        switch htype
            case {'symtensor'}
                % H = sum_{j,k} [ w_j' * h1_{j,k} * w_k ]
                % ^ where j, k are classes/modalities
                ind = spm_matcomp('SymIndices', size(h1, 4));
                w = reshape(w, [prod(dim_lattice) dim_classes dim_latent]);
                for j=1:dim_classes
                    for k=1:dim_classes
                        wk = reshape(w(:,k,:), [], dim_latent);
                        hjk = h1(:,:,:,ind(j, k));
                        % > h1_{j,k} * w_k
                        hjk = bsxfun(@times, hjk(:), wk);
                        clear wk
                        wj = reshape(w(:,j,:), [], dim_latent);
                        % > w_j' * h1_{j,k} * w_ka
                        h = h + double(wj' * hjk);
                        clear hjk wj
                    end
                end
            case {'diagonal'}
                % H = sum_k [ w_k' * h1_k * w_k ]
                % ^ where k is a class/modality
                w = reshape(w, [prod(dim_lattice) dim_classes dim_latent]);
                for k=1:dim_classes
                    wk = reshape(w(:,k,:), [], dim_latent);
                    hk = h1(:,:,:,k);
                    h = h + double(wk' * bsxfun(@times, hk(:), wk));
                    clear wk hk
                end
            otherwise
                error('Unknown hessian type.');
        end
    else
        h = [];
    end

    if do_hessian && ~do_gradient
        [g, h] = deal(h, g);
    end
    
end
