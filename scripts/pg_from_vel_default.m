function opt = pg_from_vel_default(opt)
% FORMAT opt = pg_from_vel_default(opt)
%
% Set default options for the principal geodesic model

    if nargin < 1
        opt = struct;
    end
    
    % ---------------------------------------------------------------------
    % Global parameters
    % ---------------------------------------------------------------------

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
    if ~isfield(opt, 'Az0')
        opt.Az0 = eye(opt.K);
    end
    if ~isfield(opt, 'nz0')
        opt.nz0 = max(opt.K, 10);
    end
    if ~isfield(opt, 'lam0')
        opt.lam0 = 10;
    end
    if ~isfield(opt, 'nlam0')
        opt.nl0 = 10;
    end

    
    % --- Check n0
    if 0 < opt.nz0 && opt.nz0 < opt.K
        warning(['Wishart prior: nz0 must be greater or equal to %d. ' ...
                 'Fixing it.'], opt.K)
        opt.nz0 = opt.K;
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
        'w',   'subspace.nii' ...
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
        'v',    'initial_velocity.nii',      ...
        'm',    'initial_momentum.nii'      ...
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