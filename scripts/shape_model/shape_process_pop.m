function model = shape_process_pop(dat, model, opt)
% FORMAT dat = shape_process_pop(dat, model, opt)
% dat   - Subject-specific data
% model - Model-specific data
% opt   - Options
%
% `model` and `opt` are structures that can either be in memory or on disk 
% in the form of a mat file. In the latter case, it is read and, if needed,
% written back.
% `dat` is either a structure array, a cell of structures or a cell of mat
% files. It is also read/written if needed.
%--------------------------------------------------------------------------
% Perform population-specific updates.
% This function can be run as an independent job on a cluster.
%__________________________________________________________________________
% Copyright (C) 2018 Wellcome Centre for Human Neuroimaging

    % =====================================================================
    % Read input from disk (if needed)
    [dat, ~, model, modelpath, opt] = fileToStruct(dat, model, opt);
    
    if model.converged
        if opt.verbose
            fprintf('Already converged. Unnecessary iteration.\n');
        end
        return
    end
        
    % =====================================================================
    % GENERAL TRACKING
    % > Compute LB gain (eventually, performin a moving average to smooth
    %   changes due to stochastic trace and log-determinant approximations.
    % > If LB converged, activate new components (or exit)
    
    model.emit = model.emit + 1;

    N = numel(model.lb.lb.gainlist);
    moving_gain = mean(abs(model.lb.lb.gainlist(N:-1:max(1,N-opt.lb.moving+1))));
    if moving_gain < opt.lb.threshold
        if opt.optimise.q.q && ~model.q.active
            model.q.active = true;
            fprintf('%10s | %10s\n', 'Activate', 'Affine');
        elseif opt.optimise.v.v && ~model.v.active
            model.v.active = true;
            fprintf('%10s | %10s\n', 'Activate', 'Velocity');
        elseif (opt.optimise.pg.w || opt.optimise.z.z) && ~model.pg.active
            model.pg.active = true;
            fprintf('%10s | %10s\n', 'Activate', 'PG');
        else
            fprintf('Converged :D\n');
            model.converged = true;
            model = exit_function(model, ondisk.model);
            return
        end
    end
    
    if opt.ui.verbose
        fprintf(['%10s | %10d | ' repmat('=',1,50) ' |\n'], 'VEM', model.emit);
    end

    % =====================================================================
    % AFFINE PRIOR
    if model.q.active
        % Aggregate sufficient statistics
        model = aggregateAffine(dat, model, opt);
        if opt.optimise.q.A && opt.q.Mr
            
            % Update
            model = updateAffinePrior(model, opt);
            % Lower bound
            model = lbAffinePrior(model, opt);
            if iscell(dat)
                parfor (n=1:numel(dat), opt.par.within_main)
                    dat{n} = lbAffine(dat{n}, model, opt);
                end
            else
                parfor (n=1:numel(dat), opt.par.within_main)
                    dat(n) = lbAffine(dat(n), model, opt);
                end
            end
            model = aggregateAffine(dat, model, opt);
            model = updateLowerBound(model);
        end
    end
    
    % =====================================================================
    % VELOCITY RESIDUAL PRECISION
    if model.v.active
        % Aggregate sufficient statistics
        model = aggregateVelocity(dat, model, opt);
        if opt.optimise.v.l
            
            % Update
            model = updateResidualPrecision(dat, model, opt);
            % Lower bound
            model = lbResidualPrecision(dat, model, opt);
            if iscell(dat)
                parfor (n=1:numel(dat), opt.par.within_main)
                    dat{n} = lbVelocityShape(dat{n}, model, opt);
                end
            else
                parfor (n=1:numel(dat), opt.par.within_main)
                    dat(n) = lbVelocityShape(dat(n), model, opt);
                end
            end
            model = aggregateVelocity(dat, model, opt);
            model = updateLowerBound(model);
        end
    end
    
    % =====================================================================
    % PRINCIPAL SUBSPACE
    if model.pg.active
        % Aggregate sufficient statistics
        model = aggregateLatent(dat, model, opt);
        if opt.optimise.pg.w
            
            % Update
            model = updateSubspace(dat, model, opt);
            % Lower bound
            model = lbSubspace(model, opt);
            if iscell(dat)
                parfor (n=1:numel(dat), opt.par.within_main)
                    dat{n} = lbVelocityShape(dat{n}, model, opt);
                end
            else
                parfor (n=1:numel(dat), opt.par.within_main)
                    dat(n) = lbVelocityShape(dat(n), model, opt);
                end
            end
            model = aggregateVelocity(dat, model, opt);
            model = updateLowerBound(model);
        end
    end
    
    % =====================================================================
    % ORTHOGONALISATION / LATENT PRIOR
    if model.pg.active 
        orthogonalise = opt.optimise.z.z && opt.optimise.z.A && opt.optimise.pg.w;
        [~,p] = chol(model.pg.ww); % Check pos-def
        if orthogonalise && p == 0
            
            % Update
            model = orthogonaliseSubspace(model, opt);
            % Lower bound
            model = lbLatentPrior(model, opt);
            model = lbSubspace(model, opt);
            if iscell(dat)
                parfor (n=1:numel(dat), opt.par.within_main)
                    dat{n} = lbLatent(dat{n}, model, opt);
                end
            else
                parfor (n=1:numel(dat), opt.par.within_main)
                    dat(n) = lbLatent(dat(n), model, opt);
                end
            end
            model = aggregateLatent(dat, model, opt);
            model = updateLowerBound(model);
            
        elseif opt.optimise.z.A
            
            % Update
            model = updateLatentPrior(model, opt);
            % Lower bound
            model = lbLatentPrior(model, opt);
            if iscell(dat)
                parfor (n=1:numel(dat), opt.par.within_main)
                    dat{n} = lbLatent(dat{n}, model, opt);
                end
            else
                parfor (n=1:numel(dat), opt.par.within_main)
                    dat(n) = lbLatent(dat(n), model, opt);
                end
            end
            model = aggregateLatent(dat, model, opt);
            model = updateLowerBound(model);
            
        end
    end
    
    
    % =====================================================================
    % TEMPLATE
    if opt.f.N && opt.optimise.tpl.a
        switch lower(opt.tpl.update)
            case 'map'
                model = aggregateTemplateGradHess(dat, model, opt);
                model = updateTemplate(model, opt);
            case 'ml'
                model = updateTemplateML(dat, model, opt);
        end
    end
    
    
    % =====================================================================
    % Write (if needed)
    model = structToFile(model, modelpath);
    
    % =====================================================================
    % Save everything (to allow starting from a previous state)
    if ~isempty(opt.fnames.result)
        % Ensure nifti headers are ok
        createAllNifti(dat, model, opt);
        % Write workspace
        ftrack = opt.ui.ftrack;
        opt.ui.ftrack = nan;
        save(fullfile(opt.dir.model, opt.fnames.result), ...
             'model', 'dat', 'opt');
        opt.ui.ftrack = ftrack;
    end
end