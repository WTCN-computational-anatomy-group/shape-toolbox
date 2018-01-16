function opt = pgva_model_default(opt)
% FORMAT opt = pgva_model_default(opt)
%
% Set default options for the principal geodesic model

    if nargin < 1
        opt = struct;
        opt.v.N = 0;
        opt.f.N = 0;
        opt.model.dim = 3;
        opt.model.nc  = 1;
    end
    
    % =====================================================================
    %    MODEL
    % =====================================================================
    
    % ---------------------------------------------------------------------
    % Data model
    % ---------------------------------------------------------------------
    
    if ~isfield(opt, 'model')
        opt.model = struct;
        opt.model = struct('name', 'normal', 'sigma2', 1, 'b', 1);
    end
    if ~isfield(opt.model, 'name')
        opt.model.name = 'normal';
    end
    if strcmpi(opt.model.name, 'normal') && ~isfield(opt.model, 'sigma2')
        opt.model.sigma2 = 1;
    end
    if strcmpi(opt.model.name, 'laplace') && ~isfield(opt.model, 'b')
        opt.model.b = 1;
    end
    opt.tpl.cat = any(strcmpi(opt.model.name, {'bernoulli', 'binomial', 'categorical', 'multinomial'}));
    if opt.tpl.cat
        if opt.model.nc == 1
            opt.model.name = 'bernoulli';
        else
            opt.model.name = 'categorical';
        end
    end
    
    % ---------------------------------------------------------------------
    % Subspace
    % ---------------------------------------------------------------------
    
    if ~isfield(opt, 'pg')
        opt.pg = struct;
    end
    if ~isfield(opt.pg, 'K')
        % Number of principal components
        opt.pg.K = max(32, opt.v.N + opt.f.N);
    end
    if ~isfield(opt.pg, 'prm')
        % Parameters of the differential operator
        % 1) Absolute displacement
        % 2) Membrane energy
        % 3) Bending energy
        % 4) Linear elastic energy - length change
        % 5) Linear elastic energy - volume change
        opt.pg.prm = [0.0001 0.001 0.2 0.05 0.2];
    end
    if ~isfield(opt.pg, 'bnd')
        % Boundary condition of the differential operator
        % 0 = Circulant (translation invariant)
        % 1 = Neumann   (null displacement on borders)
        % 2 = Dirichlet (null derivatives on borders)
        % 3 = Sliding   (rotation invariant)
        opt.pg.bnd = 0;
    end
    
    % ---------------------------------------------------------------------
    % Template
    % ---------------------------------------------------------------------
    
    if ~isfield(opt, 'tpl')
        opt.tpl = struct;
    end
%     error('Smooth prior not used for now')
%     if ~isfield(opt.tpl, 'prm')
%         % Parameters of the differential operator
%         % 1) Absolute displacement
%         % 2) Membrane energy
%         % 3) Bending energy
%         opt.tpl.prm = [1e-3  1e-1 0];
%     end
    if ~isfield(opt.tpl, 'bnd')
        % Boundary condition for interpolation
        % 0 = Mirror
        % 1 = Circulant (translation invariant)
        opt.tpl.bnd = 1;
    end
    if ~isfield(opt.tpl, 'itrp')
        % Interpolation order when warping the template to subjects
        opt.tpl.itrp = 1;
    end
    
    % ---------------------------------------------------------------------
    % Velocity
    % ---------------------------------------------------------------------
    
    if ~isfield(opt, 'v')
        opt.v = struct;
    end
    if ~isfield(opt.v, 'l0')
        % Prior value for lambda
        % Lambda is the anatomical noise not capture by the
        % principal subspace.
        opt.v.l0 = 17;
    end
    if ~isfield(opt.v, 'n0')
        % Prior degrees of freedom for lambda
        opt.v.n0 = 10;
    end
    
    % ---------------------------------------------------------------------
    % Latent coordnates
    % ---------------------------------------------------------------------
    
    if ~isfield(opt, 'z')
        opt.z = struct;
    end
    if ~isfield(opt.z, 'init')
        % Method used to initialise latent coordinates.
        % 'auto' = If a subspace is provided, 'zero', else, 'rand'.
        % 'zero' = Start from zero coordinates
        % 'rand' = Start from random coordinates
        opt.z.init = 'auto';
    end
    if ~isfield(opt.z, 'A0')
        % Prior value for the latent precision matrix
        opt.z.A0 = eye(opt.pg.K);
    end
    if ~isfield(opt.z, 'n0')
        % Prior DF for the latent precision matrix
        % Must be > K
        opt.z.n0 = opt.pg.K;
    end
    
    % ---------------------------------------------------------------------
    % Affine coordnates
    % ---------------------------------------------------------------------
    
    if ~isfield(opt, 'q')
        opt.q = struct;
    end
    if ~isfield(opt.q, 'B')
        % Affine transformation space
        % 'translation'/'rotation'/'rigid'/'smilitude'/'affine'
        opt.q.B = 'rigid';
    end
    if ischar(opt.q.B) || isscalar(opt.q.B)
        if opt.model.dim == 2
            dim = '2d';
        else
            dim = '3d';
        end
        [opt.q.B, opt.q.rind] = affine_basis(opt.q.B, dim);
    end
    opt.q.M  = size(opt.q.B, 3);
    opt.q.Mr = numel(opt.q.rind);
    if ~isfield(opt.q, 'A0')
        % Prior value for the affine precision matrix
        opt.q.A0 = eye(opt.q.Mr);
    end
    if ~isfield(opt.q, 'n0')
        % Prior DF for the affine precision matrix
        % Must be > Mr
        opt.q.n0 = opt.q.Mr;
    end
    if ~isfield(opt.q, 'hapx')
        % Use an approximation for the affine Hessian
        opt.q.hapx = true;
    end
    
    % =====================================================================
    %    PROCESSING
    % =====================================================================
    
    % ---------------------------------------------------------------------
    % Iterations
    % ---------------------------------------------------------------------
    
    if ~isfield(opt, 'iter')
        opt.iter = struct;
    end
    if ~isfield(opt.iter, 'em')
        % Maximum number of EM iterations
        % It should not be useful in theory, as convergence is assessed by
        % lower bound gain.
        opt.iter.em = 1000;
    end
    if ~isfield(opt.iter, 'gn')
        % Maximum number of Gauss-Newton iterations.
        % Setting it higher than 1 may help, since it improves the Laplace
        % approximation.
        opt.iter.gn = 1;
    end
    if ~isfield(opt.iter, 'ls')
        % Maximum number of line search iterations.
        % They are used to avoid GN overshooting.
        opt.iter.ls = 6;
    end
    if ~isfield(opt.iter, 'itg')
        % Number of integration steps for shooting equations
        % (NaN = try to guess something efficient)
        opt.iter.itg = nan;
    end
    
    % ---------------------------------------------------------------------
    % Lower bound
    % ---------------------------------------------------------------------
    
    if ~isfield(opt, 'lb')
        opt.lb = struct;
    end
    if ~isfield(opt.lb,  'threshold')
        % Gain threshold under which convergence is assumed
        opt.lb.threshold = 1e-5;
    end
    
    % ---------------------------------------------------------------------
    % Split processing
    % ---------------------------------------------------------------------
    
    if ~isfield(opt, 'split')
        opt.split = struct;
    end
    if ~isfield(opt.split, 'par')
        % Number of workers for parallel processing
        % 0   = no parallelisation
        % inf = Automatically set by Matlab
        opt.split.par = inf;
    end
    if ~isfield(opt.split, 'loop')
        % Along which dimension should processing be parallelised
        % 'none'/'subject'/'slice'
        opt.split.loop = 'subject';
    end
    if ~isfield(opt.split, 'batch')
        % Batch size when distributing subjects
        if opt.split.par > 0 && isfinite(opt.split.par)
            opt.split.batch = opt.split.par;
        elseif opt.split.par == 0
            opt.split.batch = 1;
        else
            myCluster = parcluster;
            opt.split.batch = myCluster.NumWorkers;
            clear myCluster
        end
    end
    
    % ---------------------------------------------------------------------
    % User interaction
    % ---------------------------------------------------------------------
    
    if ~isfield(opt, 'ui')
        opt.ui = struct;
    end
    if ~isfield(opt.ui, 'verbose')
        % Talk
        opt.ui.verbose = true;
    end
    if ~isfield(opt.ui, 'debug')
        % Talk more
        opt.ui.debug = false;
    end
    if ~isfield(opt.ui, 'ftrack')
        % Figure object where to plot the lower bound
        opt.ui.ftrack = gcf;
    end
    
    
    % ---------------------------------------------------------------------
    % Optimise/Fixed
    % ---------------------------------------------------------------------
    
    if ~isfield(opt, 'optimise')
        opt.optimise = struct;
    end
    if ~isfield(opt.optimise, 'pg')
        opt.optimise.pg = struct;
    end
    if ~isfield(opt.optimise.pg, 'w')
        opt.optimise.pg.w = true;
    end
    if ~isfield(opt.optimise, 'q')
        opt.optimise.q = struct;
    end
    if ~isfield(opt.optimise.q, 'A')
        opt.optimise.q.A = true;
    end
    if ~isfield(opt.optimise.q, 'q')
        opt.optimise.q.q = true;
    end
    if ~isfield(opt.optimise, 'z')
        opt.optimise.z = struct;
    end
    if ~isfield(opt.optimise.z, 'A')
        opt.optimise.z.A = true;
    end
    if ~isfield(opt.optimise.z, 'z')
        opt.optimise.z.z = true;
    end
    if ~isfield(opt.optimise, 'v')
        opt.optimise.v = struct;
    end
    if ~isfield(opt.optimise.v, 'l')
        opt.optimise.v.l = true;
    end
    if ~isfield(opt.optimise.v, 'v')
        opt.optimise.v.v = true;
    end
    if ~isfield(opt.optimise, 'tpl')
        opt.optimise.tpl = struct;
    end
    if ~isfield(opt.optimise.tpl, 'a')
        opt.optimise.tpl.a = true;
    end
    
    % =====================================================================
    %    CHECK SOME STUFF
    % =====================================================================
    
    % --- Check Wishart DF
    if 0 < opt.z.n0 && opt.z.n0 < opt.pg.K
        warning(['Wishart prior: z.n0 must be greater or equal to %d. ' ...
                 'Fixing it.'], opt.pg.K)
        opt.z.n0 = opt.K;
    end
    if 0 < opt.q.n0 && opt.q.n0 < opt.q.Mr
        warning(['Wishart prior: q.n0 must be greater or equal to %d. ' ...
                 'Fixing it.'], opt.q.Mr)
        opt.nq0 = opt.q.Mr;
    end
    % --- Number of PCs
    if opt.pg.K > opt.v.N + opt.f.N
        warning(['No point learnign more than %d principal components. ' ...
                 'Fixing it.'], opt.v.N + opt.f.N)
        opt.pg.K = opt.v.N + opt.f.N;
    end
        
    
    % =====================================================================
    %    I/O
    % =====================================================================
    
    % ---------------------------------------------------------------------
    % On disk / On memory
    % ---------------------------------------------------------------------
    
    if ~isfield(opt, 'ondisk')
        opt.ondisk = struct;
    end
    
    % --- Model
    
    if ~isfield(opt.ondisk, 'model')
        opt.ondisk.model = struct;
    end
    
    if ~isfield(opt.ondisk.model, 'pg')
        opt.ondisk.model.pg = struct;
    end
    if ~isfield(opt.ondisk.model.pg, 'w')
        opt.ondisk.model.pg.w = true;
    end
    
    if ~isfield(opt.ondisk.model, 'tpl')
        opt.ondisk.model.tpl = struct;
    end
    if ~isfield(opt.ondisk.model.tpl, 'mu')
        opt.ondisk.model.tpl.mu = true;
    end
    if ~isfield(opt.ondisk.model.tpl, 'gmu')
        opt.ondisk.model.tpl.gmu = true;
    end
    if ~isfield(opt.ondisk.model.tpl, 'a')
        opt.ondisk.model.tpl.a = true;
    end
        
    % --- Subjects
        
    if ~isfield(opt.ondisk, 'dat')
        opt.ondisk.dat = struct;
    end
    
    if ~isfield(opt.ondisk.dat, 'v')
        opt.ondisk.dat.v = struct;
    end
    if ~isfield(opt.ondisk.dat.v, 'v')
        opt.ondisk.dat.v.v = true;
    end
    if ~isfield(opt.ondisk.dat.v, 'm')
        opt.ondisk.dat.v.m = true;
    end
    if ~isfield(opt.ondisk.dat.v, 'iphi')
        opt.ondisk.dat.v.iphi = false;
    end
    if ~isfield(opt.ondisk.dat.v, 'ipsi')
        opt.ondisk.dat.v.ipsi = false;
    end
    if ~isfield(opt.ondisk.dat.v, 'phi')
        opt.ondisk.dat.v.phi = false;
    end
    if ~isfield(opt.ondisk.dat.v, 'jac')
        opt.ondisk.dat.v.jac = false;
    end
    if ~isfield(opt.ondisk.dat.v, 'r')
        opt.ondisk.dat.v.r = false;
    end
    if ~isfield(opt.ondisk.dat.v, 'g')
        opt.ondisk.dat.v.g = false;
    end
    if ~isfield(opt.ondisk.dat.v, 'h')
        opt.ondisk.dat.v.h = false;
    end
    
    if ~isfield(opt.ondisk.dat, 'tpl')
        opt.ondisk.dat.tpl = struct;
    end
    if ~isfield(opt.ondisk.dat.tpl, 'wa')
        opt.ondisk.dat.tpl.wa = true;
    end
    if ~isfield(opt.ondisk.dat.tpl, 'wmu')
        opt.ondisk.dat.tpl.wmu = true;
    end
    
    if ~isfield(opt.ondisk.dat, 'f')
        opt.ondisk.dat.f = struct;
    end
    if ~isfield(opt.ondisk.dat.f, 'f')
        opt.ondisk.dat.f.f = true;
    end
    if ~isfield(opt.ondisk.dat.f, 'pf')
        opt.ondisk.dat.f.pf = true;
    end
    if ~isfield(opt.ondisk.dat.f, 'c')
        opt.ondisk.dat.f.c = true;
    end
    
    % ---------------------------------------------------------------------
    % Directory
    % ---------------------------------------------------------------------
    
    if ~isfield(opt, 'dir')
        opt.dir = struct;
    end
    if ~isfield(opt.dir, 'model')
        opt.dir.model = '.';
    end
    if ~isfield(opt.dir, 'dat')
        % If empty -> next to the input
        % Else, create one directory per subject at the specified location
        opt.dir.dat = '';
    end
    
    if ~isfield(opt, 'fnames')
        opt.fnames = struct;
    end
    
    % ---------------------------------------------------------------------
    % Model
    % ---------------------------------------------------------------------
    
    if ~isfield(opt.fnames, 'result')
        opt.fnames.result = 'pgva_model.mat';
    end
    
    if ~isfield(opt.fnames, 'model')
        opt.fnames.model = struct;
    end
    
    if ~isfield(opt.fnames.model, 'pg')
        opt.fnames.model.pg = struct;
    end
    if ~isfield(opt.fnames.model.pg, 'w')
        opt.fnames.model.pg.w = 'subspace.nii';
    end
    
    if ~isfield(opt.fnames.model, 'tpl')
        opt.fnames.model.tpl = struct;
    end
    if ~isfield(opt.fnames.model.tpl, 'mu')
        opt.fnames.model.tpl.mu = 'template.nii';
    end
    if ~isfield(opt.fnames.model.tpl, 'gmu')
        opt.fnames.model.tpl.gmu = 'grad_template.nii';
    end
    if ~isfield(opt.fnames.model.tpl, 'a')
        opt.fnames.model.tpl.a = 'log_template.nii';
    end
    
    % ---------------------------------------------------------------------
    % Subjects
    % ---------------------------------------------------------------------
    
    if ~isfield(opt.fnames, 'dat')
        opt.fnames.dat = struct;
    end
    
    if ~isfield(opt.fnames.dat, 'v')
        opt.fnames.dat.v = struct;
    end
    if ~isfield(opt.fnames.dat.v, 'v')
        if isempty(opt.dir.dat)
            opt.fnames.dat.v.v = 'v_';
        else
            opt.fnames.dat.v.v = 'velocity.nii';
        end
    end
    if ~isfield(opt.fnames.dat.v, 'm')
        if isempty(opt.dir.dat)
            opt.fnames.dat.v.m = 'm_';
        else
            opt.fnames.dat.v.m = 'momentum.nii';
        end
    end
    if ~isfield(opt.fnames.dat.v, 'iphi')
        if isempty(opt.dir.dat)
            opt.fnames.dat.v.iphi = 'iphi_';
        else
            opt.fnames.dat.v.iphi = 'inverse_diffeo.nii';
        end
    end
    if ~isfield(opt.fnames.dat.v, 'ipsi')
        if isempty(opt.dir.dat)
            opt.fnames.dat.v.ipsi = 'ipsi_';
        else
            opt.fnames.dat.v.ipsi = 'inverse_transform.nii';
        end
    end
    if ~isfield(opt.fnames.dat.v, 'phi')
        if isempty(opt.dir.dat)
            opt.fnames.dat.v.phi = 'phi_';
        else
            opt.fnames.dat.v.phi = 'forward_diffeo.nii';
        end
    end
    if ~isfield(opt.fnames.dat.v, 'jac')
        if isempty(opt.dir.dat)
            opt.fnames.dat.v.jac = 'j_';
        else
            opt.fnames.dat.v.jac = 'jacobian_diffeo.nii';
        end
    end
    if ~isfield(opt.fnames.dat.v, 'r')
        if isempty(opt.dir.dat)
            opt.fnames.dat.v.r = 'r_';
        else
            opt.fnames.dat.v.r = 'residual_field.nii';
        end
    end
    if ~isfield(opt.fnames.dat.v, 'g')
        if isempty(opt.dir.dat)
            opt.fnames.dat.v.g = 'gv_';
        else
            opt.fnames.dat.v.g = 'gradient_velocity.nii';
        end
    end
    if ~isfield(opt.fnames.dat.v, 'h')
        if isempty(opt.dir.dat)
            opt.fnames.dat.v.h = 'hv_';
        else
            opt.fnames.dat.v.h = 'hessian_velocity.nii';
        end
    end
    
    if ~isfield(opt.fnames.dat, 'tpl')
        opt.fnames.dat.tpl = struct;
    end
    if ~isfield(opt.fnames.dat.tpl, 'wa')
        if isempty(opt.dir.dat)
            opt.fnames.dat.tpl.wa = 'wa_';
        else
            opt.fnames.dat.tpl.wa = 'warped_log_template.nii';
        end
    end
    if ~isfield(opt.fnames.dat.tpl, 'wmu')
        if isempty(opt.dir.dat)
            opt.fnames.dat.tpl.wmu = 'wmu_';
        else
            opt.fnames.dat.tpl.wmu = 'warped_template.nii';
        end
    end
    
    if ~isfield(opt.fnames.dat, 'f')
        opt.fnames.dat.f = struct;
    end
    if ~isfield(opt.fnames.dat.f, 'f')
        if isempty(opt.dir.dat)
            opt.fnames.dat.f.f = 'f_';
        else
            opt.fnames.dat.f.f = 'observed_image.nii';
        end
    end
    if ~isfield(opt.fnames.dat.f, 'pf')
        if isempty(opt.dir.dat)
            opt.fnames.dat.f.pf = 'pf_';
        else
            opt.fnames.dat.f.pf = 'pushed_image.nii';
        end
    end
    if ~isfield(opt.fnames.dat.f, 'c')
        if isempty(opt.dir.dat)
            opt.fnames.dat.f.c = 'c_';
        else
            opt.fnames.dat.f.c = 'pushed_count.nii';
        end
    end
    
end