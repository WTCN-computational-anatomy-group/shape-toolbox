function opt = pgr_registration(opt)
% FORMAT opt = pgr_registration(opt)
%
% Performs diffeomorphic registration of a template towards a target image,
% based on a "principal geodesic" model + a residual field accounting for
% the "geodesic noise".
%
% The input option structure (opt) should at least contain the fields:
% fnames.f  / dat.f  - the target image (as a file or array)
% fnames.mu / dat.mu - the template image (as a file or array)
% fnames.w  / dat.w  - the principal sibspace (as a file or array)
%
% The following parameters can be overriden by specifying them in the
% input option structure:
% model        - Generative model [struct('name', 'normal')]
% sigma        - Variance captured by the residual field [1]
% itrp         - Interpolation order [1]
% bnd          - Boundary conditions [1]
% gnit         - Number of Gauss-Newton iterations [1]
% lsit         - Number of line search iterations [6]
% regz         - Precision matrix for the latent coordinates [W'LW]
% itgr         - Number of integration steps for geodesic shooting [auto]
% prm          - Parameters of the geodesic differential operator
%                [0.0001 0.001 0.2 0.05 0.2]
% verbose      - Talk during line search [true]
% debug        - Further debuging talk [false]
% loop         - How to split array processing [auto]
% par          - Parallelise processing [true/auto]
% directory    - Directory where to store result arrays ['.']
% fnames       - Structure of filenames for all temporary arrays:
%                (f, mu, w, gmu, r, v, ipsi, iphi, pf, c, (a))
%
% The output option structure contains the additional field:
% dat       - Structure of file_arrays and arrays storing all results
%             (z, r, v, llm, llq, f, mu, gmu, iphi, iphi, pf, c, (a))
% It can also be used as input to specify initial starting points.

    % ---------------------------------------------------------------------
    %    Default parameters
    % ---------------------------------------------------------------------

    % Global parameters
    % -----------------
    if ~isfield(opt, 'model')
        opt.model = struct('name', 'normal', 'sigma2', 1);
    end
    opt.tpm = any(strcmpi(opt.model.name, {'bernoulli', 'binomial', 'categorical', 'multinomial'}));
    if ~isfield(opt, 'sigma')
        opt.sigma = 1;
    end
    if ~isfield(opt, 'itrp')
        opt.itrp = 1;
    end
    if ~isfield(opt, 'bnd')
        opt.bnd = 1;
    end
    if ~isfield(opt, 'gnit')
        opt.gnit = 1;
    end
    if ~isfield(opt, 'lsit')
        opt.lsit = 6;
    end
    if ~isfield(opt, 'regz')
        opt.regz = [];
    end
    if ~isfield(opt, 'itgr')
        opt.itgr = nan;
    end
    if ~isfield(opt, 'prm')
        opt.prm = [0.0001 0.001 0.2 0.05 0.2];
    end
    if ~isfield(opt, 'debug')
        opt.debug = false;
    end
    if ~isfield(opt, 'verbose')
        opt.verbose = true;
    end
    if ~isfield(opt, 'par')
        opt.par = true;
    end
    if ~isfield(opt, 'loop')
        opt.loop = '';
    end
    
    % Files
    % -----
    if ~isfield(opt, 'directory')
        opt.directory = '.';
    end
    if ~isfield(opt, 'fnames')
        opt.fnames = struct;
    end
    if ~isfield(opt.fnames, 'f')
        opt.fnames.f = fullfile(opt.directory, 'testimage.nii');
    end
    if ~isfield(opt.fnames, 'mu')
        opt.fnames.mu = fullfile(opt.directory, 'template.nii');
    end
    if ~isfield(opt.fnames, 'w')
        opt.fnames.w = fullfile(opt.directory, 'subspace.nii');
    end
    if ~isfield(opt.fnames, 'wmu')
        opt.fnames.wmu = fullfile(opt.directory, 'warped_template.nii');
    end
    if ~isfield(opt.fnames, 'gmu')
        opt.fnames.gmu = fullfile(opt.directory, 'grad_template.nii');
    end
    if ~isfield(opt.fnames, 'iphi')
        opt.fnames.iphi = fullfile(opt.directory, 'inverse_diffeo_warp.nii');
    end
    if ~isfield(opt.fnames, 'ipsi')
        opt.fnames.ipsi = fullfile(opt.directory, 'inverse_complete_warp.nii');
    end
    if ~isfield(opt.fnames, 'r')
        opt.fnames.r = fullfile(opt.directory, 'residual_field.nii');
    end
    if ~isfield(opt.fnames, 'v')
        opt.fnames.v = fullfile(opt.directory, 'initial_velocity.nii');
    end
    if ~isfield(opt.fnames, 'pf')
        opt.fnames.pf = fullfile(opt.directory, 'pushed.nii');
    end
    if ~isfield(opt.fnames, 'c')
        opt.fnames.c = fullfile(opt.directory, 'count.nii');
    end
    if opt.tpm
        if ~isfield(opt.fnames, 'a')
            opt.fnames.a = fullfile(opt.directory, 'logtemplate.nii');
        end
    end
    
    % Initial state
    % -------------
    if ~isfield(opt, 'dat')
        opt.dat = struct;
    end
    
    % -- Target image
    if ~isfield(opt.dat, 'f')
        if ~exist(opt.fnames.f, 'file')
            error('A target image must be provided')
        end
        n = nifti(opt.fnames.f);
        opt.dat.f = n.dat;
    end
    if ~isfield(opt.dat, 'Mf')
        if isa(opt.dat.f, 'file_array') && endsWith(opt.dat.f.fname, '.nii')
            n = nifti(opt.dat.f.fname);
            opt.dat.Mf = n.mat0;
        else
            opt.dat.Mf = eye(4);
        end
    end
    % -- Template: Tissue probability maps case
    if opt.tpm
        if ~isfield(opt.dat, 'a')
            if ~exist(opt.fnames.a, 'file')
                error('A log-template image must be provided')
            end
            n = nifti(opt.fnames.a);
            opt.dat.a = n.dat;
        end
        if ~isfield(opt.dat, 'Mmu')
            if isa(opt.dat.a, 'file_array') && endsWith(opt.dat.a.fname, '.nii')
                n = nifti(opt.dat.a.fname);
                opt.dat.Mmu = n.mat0;
            else
                opt.dat.Mmu = eye(4);
            end
        end
        if ~isfield(opt.dat, 'mu')
            opt.dat.mu = file_array(opt.fnames.mu, size(opt.dat.a), 'float32');
            n = nifti;
            n.dat = opt.dat.mu;
            n.mat0 = opt.dat.Mmu;
            create(n);
        end
    % -- Template: Mean intensity case
    else
        if ~isfield(opt.dat, 'mu')
            if ~exist(opt.fnames.mu, 'file')
                error('A template image must be provided')
            end
            n = nifti(opt.fnames.mu);
            opt.dat.mu = n.dat;
        end
        if ~isfield(opt.dat, 'Mmu')
            if isa(opt.dat.mu, 'file_array') && endsWith(opt.dat.mu.fname, '.nii')
                n = nifti(opt.dat.mu.fname);
                opt.dat.Mmu = n.mat0;
            else
                opt.dat.Mmu = eye(4);
            end
        end
    end
    % -- Principal subspace
    if ~isfield(opt.dat, 'w')
        if ~exist(opt.fnames.w, 'file')
            error('A subspace must be provided')
        end
        n = nifti(opt.fnames.w);
        opt.dat.w = n.dat;
    end
    if ~isfield(opt.dat, 'z')
        opt.dat.z = zeros(size(opt.dat.w, 5), 1, 'double');
    end
    % -- Gradient template
    if ~isfield(opt.dat, 'gmu')
        dim = [size(opt.dat.mu) 1 1];
        dim = dim(1:4);
        opt.dat.gmu = file_array(opt.fnames.gmu, [dim 3], 'float32');
        n = nifti;
        n.dat = opt.dat.gmu;
        n.mat0 = opt.dat.Mmu;
        create(n);
    end
    % -- Warped template
    if ~isfield(opt.dat, 'wmu')
        dim = [size(opt.dat.f) 1 1];
        dim = dim(1:4);
        opt.dat.wmu = file_array(opt.fnames.wmu, [dim size(opt.dat.mu, 4)], 'float32');
        n = nifti;
        n.dat = opt.dat.wmu;
        n.mat0 = opt.dat.Mf;
        create(n);
    end
    % -- Transform & Pushed image
    if ~isfield(opt.dat, 'pf')
        lat = [size(opt.dat.mu) 1];
        lat = lat(1:3);
        opt.dat.pf = file_array(opt.fnames.pf, [lat size(opt.dat.f, 4)], 'float32');
        n = nifti;
        n.dat = opt.dat.pf;
        n.mat0 = opt.dat.Mmu;
        create(n);
    end
    if ~isfield(opt.dat, 'c')
        lat = [size(opt.dat.mu) 1];
        lat = lat(1:3);
        opt.dat.c = file_array(opt.fnames.c, lat, 'float32');
        n = nifti;
        n.dat = opt.dat.c;
        n.mat0 = opt.dat.Mmu;
        create(n);
    end
    if ~isfield(opt.dat, 'ipsi')
        lat = [size(opt.dat.f) 1];
        lat = lat(1:3);
        opt.dat.ipsi = file_array(opt.fnames.ipsi, [lat 3], 'float32');
        n = nifti;
        n.dat = opt.dat.ipsi;
        create(n);
    end
    if ~isfield(opt.dat, 'iphi')
        lat = [size(opt.dat.mu) 1];
        lat = lat(1:3);
        opt.dat.iphi = file_array(opt.fnames.iphi, [lat 3], 'float32');
        n = nifti;
        n.dat = opt.dat.iphi;
        create(n);
    end
    if ~isfield(opt.dat, 'v')
        lat = [size(opt.dat.mu) 1];
        lat = lat(1:3);
        opt.dat.v = file_array(opt.fnames.v, [lat 3], 'float32');
        n = nifti;
        n.dat = opt.dat.v;
        create(n);
    end
    if ~isfield(opt.dat, 'r')
        lat = [size(opt.dat.mu) 1];
        lat = lat(1:3);
        opt.dat.r = file_array(opt.fnames.r, [lat 3], 'float32');
        n = nifti;
        n.dat = opt.dat.r;
        create(n);
        opt.dat.r(:) = 0;
    end
    
    % ---------------------------------------------------------------------
    %    Initialise all variables
    % ---------------------------------------------------------------------
    
    % Exponentiate initial velocity
    % -----------------------------
    opt.dat.v = reconstructVelocity('latent', opt.dat.z, 'subspace', opt.dat.w, ...
        'residual', opt.dat.r, 'sigma', opt.sigma, ...
        'debug', opt.debug, 'output', opt.dat.v);
    opt.dat.iphi = exponentiateVelocity(opt.dat.v, 'iphi', ...
        'itgr', opt.itgr, 'vs', sqrt(sum(opt.dat.Mmu(1:3,1:3).^2)), ...
        'prm', opt.prm, 'debug', opt.debug, 'output', opt.dat.iphi);
    latf = [size(opt.dat.f) 1];
    latf= latf(1:3);
    opt.dat.ipsi = reconstructIPsi(eye(4), opt.dat.iphi, ...
        'lat', latf, 'Mf', opt.dat.Mf, 'Mmu', opt.dat.Mmu, ...
        'output', opt.dat.ipsi, 'debug', opt.debug);
    
    
    % Compute template spatial gradients (+ Build TPMs if needed)
    % -----------------------------------------------------------
    if opt.tpm
        opt.dat.gmu = templateGrad(opt.dat.a, opt.itrp, opt.bnd, ...
            'debug', opt.debug, 'output', opt.dat.gmu);
        opt.dat.mu = reconstructProbaTemplate(opt.dat.a, ...
            'loop', opt.loop, 'par', opt.par, 'debug', opt.debug, ...
            'output', opt.dat.mu);
    else
        opt.dat.gmu = templateGrad(opt.dat.mu, opt.itrp, opt.bnd, ...
            'debug', opt.debug, 'output', opt.dat.gmu);
    end
    
    % Push image to template space
    % ----------------------------
    latmu = [size(opt.dat.mu) 1];
    latmu= latmu(1:3);
    [opt.dat.pf, opt.dat.c] = pushImage(opt.dat.ipsi, opt.dat.f, latmu, ...
        'loop', opt.loop, 'par', opt.par, ...
        'output', {opt.dat.pf, opt.dat.c}, 'debug', opt.debug);
    
    % Compute initial log-likelihood
    % ------------------------------
    opt.dat.llm = llMatching(opt.model, opt.dat.mu, opt.dat.pf, opt.dat.c, ...
        'loop', opt.loop, 'par', opt.par, 'debug', opt.debug);
    
    if checkarray(opt.regz)
        opt.regz = precisionZ(opt.dat.w, sqrt(sum(opt.dat.Mmu(1:3,1:3).^2)), ...
            opt.prm, 'debug', opt.debug);
    end
    opt.dat.llz = llPriorAffine(opt.dat.z, opt.regz, 'debug', opt.debug);
    
    % ---------------------------------------------------------------------
    %    Processing
    % ---------------------------------------------------------------------

    % Gauss-Newton iterations
    % -----------------------
    for gnit = 1:opt.gnit
        
        % -----------------------------
        % 1. Optimize latent coordinate
        % -----------------------------
        
        % Compute gradient/hessian
        % ------------------------
        
        [opt.dat.gz, opt.dat.hz] = ghMatchingLatent(opt.model, ...
            opt.dat.mu, opt.dat.pf, opt.dat.c, ...
            opt.dat.gmu, opt.dat.w, ...
            'loop', opt.loop, 'par', opt.par, 'debug', opt.debug);
        
        [gz, hz] = ghPriorLatent(opt.dat.z, opt.regz, 'debug', opt.debug);
        opt.dat.gz = opt.dat.gz + gz;
        opt.dat.hz = opt.dat.hz + hz;
        clear gz hz
        
        opt.dat.hz = loadDiag(opt.dat.hz); % Additional regularisation for robustness
        
        % Compute ascent direction
        % ------------------------
        opt.dat.dz = -opt.dat.hz \ opt.dat.gz;
        
        % Line search
        % -----------
        [okz, z, llm, llz, v, iphi, pf, c, ipsi] = lsLatent(...
            opt.model, opt.dat.dz, opt.dat.z, opt.dat.v, opt.dat.llm, ...
            opt.dat.w, opt.dat.mu, opt.dat.f, ...
            'regz', opt.regz, ...
            'Mf', opt.dat.Mf, 'Mmu', opt.dat.Mmu, 'nit', opt.lsit, ...
            'itgr', opt.itgr, 'prm', opt.prm, ...
            'par', opt.par, 'verbose', opt.verbose, 'debug', opt.debug);
        
        % Store better values
        % -------------------
        if okz
            opt.dat.z       = z;
            opt.dat.llm     = llm;
            opt.dat.llz     = llz;
            opt.dat.v(:)    = v(:);
            opt.dat.iphi(:) = iphi(:);
            opt.dat.pf(:)   = pf(:);
            opt.dat.c(:)    = c(:);
            opt.dat.ipsi(:) = ipsi(:);
        end
        
        % --------------------------
        % 2. Optimize residual field
        % --------------------------
        
        % Compute gradient/hessian
        % ------------------------
        
        [opt.dat.gr, opt.dat.hr] = ghMatchingVel(opt.model, ...
            opt.dat.mu, opt.dat.pf, opt.dat.c, opt.dat.gmu, ...
            'loop', opt.loop, 'par', opt.par, 'debug', opt.debug);
        
        % Compute ascent direction
        % ------------------------
        opt.dat.dr = spm_diffeo('fmg', single(opt.dat.hr), single(opt.dat.gr), ...
            [sqrt(sum(opt.dat.Mmu(1:3,1:3).^2)) opt.prm 2 2]);
        
        % Line search
        % -----------
        [okr, r, llm, llr, iphi, pf, c, ipsi, v] = lsVelocity(...
            opt.model, opt.dat.dr, opt.dat.r, opt.dat.llm, ...
            opt.dat.mu, opt.dat.f, 'v0', opt.dat.v, ...
            'Mf', opt.dat.Mf, 'Mmu', opt.dat.Mmu, 'nit', opt.lsit, ...
            'itgr', opt.itgr, 'prm', opt.prm, ...
            'par', opt.par, 'verbose', opt.verbose, 'debug', opt.debug);
        
        % Store better values
        % -------------------
        if okr
            opt.dat.r(:)    = r(:);
            opt.dat.llm     = llm;
            opt.dat.llr     = llr;
            opt.dat.iphi(:) = iphi(:);
            opt.dat.pf(:)   = pf(:);
            opt.dat.c(:)    = c(:);
            opt.dat.ipsi(:) = ipsi(:);
            opt.dat.v(:)    = v(:);
        end
        
        
        % -----------------------------
        % 3. In no improvement, get out
        % -----------------------------
        
        if ~okz && ~okr
            break;
        end
        
    end
    
    % Update hessian for Laplace aproximation
    % ---------------------------------------
    if okz || okr
        [opt.dat.gz, opt.dat.hz] = ghMatchingLatent(opt.model, ...
            opt.dat.mu, opt.dat.pf, opt.dat.c, ...
            opt.dat.gmu, opt.dat.w, ...
            'loop', opt.loop, 'par', opt.par, 'debug', opt.debug);
        
        [gz, hz] = ghPriorLatent(opt.dat.z, opt.regz, 'debug', opt.debug);
        opt.dat.gz = opt.dat.gz + gz;
        opt.dat.hz = opt.dat.hz + hz;
        clear gz hz
    end
    
    
    % Warp template to image
    % ----------------------
    opt.dat.wmu = warp(opt.dat.ipsi, opt.dat.mu, opt.itrp, opt.bnd, ...
        'par', opt.par, 'output', opt.dat.wmu, 'debug', opt.debug);
end