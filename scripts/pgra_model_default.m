function opt = pgra_model_default(opt)
% FORMAT opt = pg_model_default(opt)
%
% Set default options for the principal geodesic model

    if nargin < 1
        opt = struct;
    end
    
    % ---------------------------------------------------------------------
    % Global parameters
    % ---------------------------------------------------------------------
    
    if ~isfield(opt, 'model')
        opt.model = struct('name', 'normal', 'sigma2', 1);
    end
    opt.tpm = any(strcmpi(opt.model.name, {'bernoulli', 'binomial', 'categorical', 'multinomial'}));
    if ~isfield(opt, 'K')
        opt.K = 32;
    end
    if ~isfield(opt, 'itrp')
        opt.itrp = 1;
    end
    if ~isfield(opt, 'bnd')
        opt.bnd = 1;
    end
    if ~isfield(opt, 'emit')
        opt.emit = 100;
    end
    if ~isfield(opt, 'lsit')
        opt.lsit = 6;
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
        opt.par = inf;
    end
    if ~isfield(opt, 'loop')
        opt.loop = 'subject';
    end
    if ~isfield(opt, 'batch')
        if opt.par > 0 && isfinite(opt.par)
            opt.batch = opt.par;
        elseif opt.par == 0
            opt.batch = 1;
        else
            myCluster = parcluster;
            opt.batch = myCluster.NumWorkers;
            clear myCluster
        end
    end
    if ~isfield(opt, 'happrox')
        opt.happrox = true;
    end
    if ~isfield(opt ,'affine_basis')
        opt.affine_basis = affine_basis(12);
    end
    if ~isfield(opt ,'affine_rind')
        opt.affine_rind = [];
    end
    if ~isfield(opt, 'nq0')
        opt.nq0 = numel(opt.affine_rind);
    end
    if ~isfield(opt, 'nz0')
        opt.nz0 = opt.K;
    end
    if ~isfield(opt, 'nlam0')
        opt.nlam0 = 10;
    end
    if ~isfield(opt, 'lambda0')
        opt.lambda0 = 10;
    end
    if ~isfield(opt, 'wpz0')
        opt.wpz0 = [1 1];
    end
    if ~isfield(opt, 'wpz')
        opt.wpz = [1 1];
    end
    if ~isfield(opt, 'armijo')
        opt.armijo = 1;
    end
    if ~isfield(opt, 'fwhm')
        opt.fwhm = 0;
    end
    
    % --- Check n0
    if 0 < opt.nz0 && opt.nz0 < opt.K
        warning(['Wishart prior: nz0 must be greater or equal to %d. ' ...
                 'Fixing it.'], opt.K)
        opt.nz0 = opt.K;
    end
    if 0 < opt.nq0 && opt.nq0 < size(opt.affine_basis, 3)
        warning(['Wishart prior: nq0 must be greater or equal to %d. ' ...
                 'Fixing it.'], size(opt.affine_basis, 3))
        opt.nq0 = size(opt.affine_basis, 3);
    end
    
    % ---------------------------------------------------------------------
    % Files
    % ---------------------------------------------------------------------
    
    % On Disk
    % -------
    
    if ~isfield(opt, 'ondisk')
        opt.ondisk = struct;
    end
    if ~isfield(opt.ondisk, 'model')
        opt.ondisk.model = struct;
    end
    if ~isfield(opt.ondisk, 'dat')
        opt.ondisk.dat = struct;
    end

    dftondisk = default_ondisk;
    
    varmodel = fieldnames(dftondisk.model);
    for i=1:numel(varmodel)
        var = varmodel{i};
        if ~isfield(opt.ondisk.model, var)
            opt.ondisk.model.(var) = dftondisk.model.(var);
        end
    end
    
    varmodel = fieldnames(dftondisk.dat);
    for i=1:numel(varmodel)
        var = varmodel{i};
        if ~isfield(opt.ondisk.dat, var)
            opt.ondisk.dat.(var) = dftondisk.dat.(var);
        end
    end
    
    % Filenames
    % ---------
    
    % --- Directory
    
    if ~isfield(opt, 'directory')
        opt.directory = '.';
    end
    if ~isfield(opt, 'fnames')
        opt.fnames = struct;
    end
    if ~isfield(opt.fnames, 'result')
        opt.fnames.result = 'pg_result.mat';
    end
    if ~isfield(opt.fnames, 'model')
        opt.fnames.model = struct;
    end
    if ~isfield(opt.fnames, 'dat')
        opt.fnames.dat = struct;
    end
    
    % --- Model
    
    dftfnames_model = struct(...
        'mu',  'template.nii',           ...
        'gmu', 'grad_template.nii',      ...
        'a',   'logtemplate.nii',        ...
        'w',   'subspace.nii',           ...
        'dw',  'subspace_ls.nii',        ...
        'gw',  'subspace_grad.nii',      ...
        'hw',  'subspace_hess.nii',      ...
        'ww',  'precision_z_smooth.nii', ...
        'A',   'precision_z.nii',        ...
        'z',   'latent_stat_z.nii',      ...
        'zz',  'latent_stat_zz.nii',     ...
        'S',   'latent_stat_cov.nii'     ...
    );

    varmodel = fieldnames(dftfnames_model);
    for i=1:numel(varmodel)
        var = varmodel{i};
        if ~isfield(opt.fnames.model, var)
            opt.fnames.model.(var) ...
                = fullfile(opt.directory, dftfnames_model.(var));
        end
    end
    
    % --- Subjects
    
    dftfnames_subjects = struct(...
        'wmu',  'warped_template.nii',       ...
        'iphi', 'inverse_diffeo_warp.nii',   ...
        'ipsi', 'inverse_complete_warp.nii', ...
        'v',    'initial_velocity.nii',      ...
        'pf',   'pushed.nii',                ...
        'c',    'count.nii',                 ...
        'gv',   'gradient_vel.nii',          ...
        'hv',   'hessian_vel.nii',           ...
        'r',    'residual_field.nii',        ...
        'z',    'latent_coordinates.nii',    ...
        'zz',   'latent_stat_zz.nii',        ...
        'S',    'latent_cov.nii'             ...
    );
    
    varmodel = fieldnames(dftfnames_subjects);
    for i=1:numel(varmodel)
        var = varmodel{i};
        if ~isfield(opt.fnames, var)
            opt.fnames.dat.(var) = cell(1, opt.N);
            for j=1:opt.N
                opt.fnames.dat.(var){j} ...
                    = fullfile(opt.directory, num2str(j), dftfnames_subjects.(var));
            end
        end
    end
    
end