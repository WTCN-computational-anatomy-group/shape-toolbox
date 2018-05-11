function [dat, model] = shape_process(dat, model, opt)
% FORMAT dat = shape_process(dat, model, opt)
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
% Perform all updates the "right" (but less efficient) way.
% It allows to accurately track the lower bound all along.
%__________________________________________________________________________
% Copyright (C) 2018 Wellcome Centre for Human Neuroimaging

    % =====================================================================
    % Read input from disk (if needed)
    [dat, datpath, model, modelpath, opt] = fileToStruct(dat, model, opt);
    
    if model.converged
        if opt.verbose
            fprintf('Already converged. Unnecessary iteration.\n');
        end
        return
    end
        
    % =====================================================================
    %   GENERAL TRACKING
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
            model = structToFile(model, modelpath);
            return
        end
    end
    
    if opt.ui.verbose, shape_ui('EM',model.emit); end

    % =====================================================================
    %   AFFINE
    if model.q.active
        % -----------------------------------------------------------------
        % Update rigid/affine transform
        if opt.optimise.q.q
            if opt.ui.verbose
                shape_ui('Title', 'Update Affine', false);
            end
            okpre = sum(toArray(dat, '.q.ok') >= 0);
            if ~isfield(opt.par.subjects.job, 'mem_affine')
                opt.par.subjects.job.mem_affine = opt.par.subjects.job.mem;
            end
            opt.par.subjects.job.mem = opt.par.subjects.job.mem_affine;
            [opt.par.subjects, dat] = distribute(opt.par.subjects, ...
                'updateAffine', 'inplace', dat, model, opt);
            opt.par.subjects.job.mem_affine = opt.par.subjects.job.mem;
            if opt.ui.verbose
                shape_ui('Title', '', false);
                okpost = sum(toArray(dat, '.q.ok') > 0);
                fprintf('success: %d / %d\n', okpost, okpre);
            end
            [~, dat] = distribute([], 'lbAffine', 'inplace', dat, model, opt);
            model = aggregateAffine(dat, model, opt);
            model = aggregateMatching(dat, model, opt);
            model = updateLowerBound(model, opt);
            shape_plot_all(model,opt);
        end
    
        % -----------------------------------------------------------------
        % Update affine prior
        if opt.optimise.q.A && opt.q.Mr
            if opt.ui.verbose
                t0 = shape_ui('Title', 'Update Affine prec', false, true);
            end
            model = updateAffinePrior(model, opt);
            model = lbAffinePrior(model, opt);
            [~, dat] = distribute([], 'lbAffine', 'inplace', dat, model, opt);
            model = aggregateAffine(dat, model, opt);
            if opt.ui.verbose
                shape_ui('PostTitle', toc(t0));
            end
            model = updateLowerBound(model, opt);
            shape_plot_all(model,opt);
        end
    end
    
    % =====================================================================
    %   VELOCITY
    if model.v.active
        % -----------------------------------------------------------------
        % Update velocity
        if opt.optimise.v.v
            if opt.ui.verbose
                shape_ui('Title', 'Update Velocity', false);
            end
            okpre = sum(toArray(dat, '.v.ok') >= 0);
            if ~isfield(opt.par.subjects.job, 'mem_velocity')
                opt.par.subjects.job.mem_velocity = opt.par.subjects.job.mem;
            end
            opt.par.subjects.job.mem = opt.par.subjects.job.mem_velocity;
            [opt.par.subjects, dat] = distribute(opt.par.subjects, ...
                'updateVelocityShape', 'inplace', dat, model, opt);
            opt.par.subjects.job.mem_velocity = opt.par.subjects.job.mem;
            if opt.ui.verbose
                shape_ui('Title', '', false);
                okpost = sum(toArray(dat, '.v.ok') > 0);
                fprintf('success: %d / %d\n', okpost, okpre);
            end
            [~, dat] = distribute([], 'lbVelocityShape', 'inplace', dat, model, opt);
            model = aggregateVelocity(dat, model, opt);
            model = aggregateMatching(dat, model, opt);
            model = updateLowerBound(model, opt);
            shape_plot_all(model,opt);
        end
        
        % =================================================================
        % Update velocity residual precision
        if opt.optimise.v.l
            if opt.ui.verbose
                t0 = shape_ui('Title', 'Update Residual prec', false, true);
                l0 = model.v.l;
            end
            model = updateResidualPrecision(model, opt);
            if opt.ui.verbose
                shape_ui('PostTitle', toc(t0));
                shape_ui('Title', '', false);
                fprintf('%5.3f -> %5.3f\n', l0, model.v.l);
            end
            model = lbResidualPrior(model, opt);
            [~, dat] = distribute([], 'lbVelocityShape', 'inplace', dat, model, opt);
            model = aggregateVelocity(dat, model, opt);
            model = updateLowerBound(model, opt);
            shape_plot_all(model,opt);
        end
    end
    
    % =====================================================================
    %   SHAPE
    if model.pg.active
        
        % -----------------------------------------------------------------
        % Update principal subspace
        if opt.optimise.pg.w
            if opt.ui.verbose
                t0 = shape_ui('Title', 'Update Subspace', false, true);
            end
            model = updateSubspace(dat, model, opt);
            model = lbSubspace(model, opt);
            [~, dat] = distribute([], 'updateResidualLB', 'inplace', dat, model, opt);
            [~, dat] = distribute([], 'lbLatent', 'inplace', dat, model, opt);
            [~, dat] = distribute([], 'lbVelocityShape', 'inplace', dat, model, opt);
            model = aggregateLatent(dat, model, opt);
            model = aggregateVelocity(dat, model, opt);
            if opt.ui.verbose
                shape_ui('PostTitle', toc(t0));
            end
            model = updateLowerBound(model, opt);
            shape_plot_all(model,opt);
        end
        
        % -----------------------------------------------------------------
        % Update latent coordinates
        if opt.optimise.z.z
            if opt.ui.verbose
                shape_ui('Title', 'Update Latent', false);
            end
            if ~isfield(opt.par.subjects.job, 'mem_latent')
                opt.par.subjects.job.mem_latent = opt.par.subjects.job.mem;
            end
            opt.par.subjects.job.mem = opt.par.subjects.job.mem_latent;
            [opt.par.subjects, dat] = distribute(opt.par.subjects, ...
                'updateLatent', 'inplace', dat, model, opt);
            opt.par.subjects.job.mem_latent = opt.par.subjects.job.mem;
            [~, dat] = distribute([], 'updateResidualLB', 'inplace', dat, model, opt);
            [~, dat] = distribute([], 'lbLatent', 'inplace', dat, model, opt);
            [~, dat] = distribute([], 'lbVelocityShape', 'inplace', dat, model, opt);
            model = aggregateLatent(dat, model, opt);
            model = aggregateVelocity(dat, model, opt);
            model = updateLowerBound(model, opt);
            shape_plot_all(model,opt);
        end
        
        % -----------------------------------------------------------------
        % Orthogonalisation / Update latent prior
        orthogonalise = opt.optimise.z.z && opt.optimise.z.A && opt.optimise.pg.w;
        [~,p] = chol(model.pg.ww);
        if orthogonalise && p == 0
            if opt.ui.verbose
                t0 = shape_ui('Title', 'Orthogonalise Subspace', false, true);
            end
            [model,Q] = orthogonaliseSubspace(model, opt);
            [~, dat] = distribute([], 'rotateLatent', Q, 'inplace', dat);
            [~, dat] = distribute([], 'lbLatent', 'inplace', dat, model, opt);
            model = aggregateLatent(dat,model,opt);
            model = lbSubspace(model, opt);
            model = lbLatentPrior(model, opt);
            if opt.ui.verbose
                shape_ui('PostTitle', toc(t0));
            end
            model = updateLowerBound(model, opt);
            shape_plot_all(model,opt);
        elseif opt.optimise.z.A
            if opt.ui.verbose
                t0 = shape_ui('Title', 'Update Latent prec', false, true);
            end
            model = updateLatentPrior(model, opt);
            model = lbLatentPrior(model, opt);
            [~, dat] = distribute([], 'lbLatent', 'inplace', dat, model, opt);
            model = aggregateLatent(dat, model, opt);
            if opt.ui.verbose
                shape_ui('PostTitle', toc(t0));
            end
            model = updateLowerBound(model, opt);
            shape_plot_all(model,opt);
        end
    end
    
    % =====================================================================
    %   TEMPLATE
    if opt.f.N && opt.optimise.tpl.a
        % -----------------------------------------------------------------
        % Compute template gradient
        if opt.ui.verbose
            shape_ui('Title', 'Update Template Grad', false);
        end
        if ~isfield(opt.par.subjects.job, 'mem_template')
            opt.par.subjects.job.mem_template = opt.par.subjects.job.mem;
        end
        opt.par.subjects.job.mem = opt.par.subjects.job.mem_template;
        [opt.par.subjects, dat] = distribute(opt.par.subjects, ...
            'updateTemplateGradHess', 'inplace', dat, model, opt);
        opt.par.subjects.job.mem_template = opt.par.subjects.job.mem;
            
        % -----------------------------------------------------------------
        % Update template
        if opt.ui.verbose
            t0 = shape_ui('Title', 'Update Template', false, true);
        end
        switch lower(opt.tpl.update)
            case 'map'
                model = aggregateTemplateGradHess(dat, model, opt);
                model = updateTemplate(model, opt);
                model = lbTemplate(model, opt);
            case 'ml'
                model = updateTemplateML(dat, model, opt);
        end
        model = updateTemplateDerivatives(model, opt);
        if opt.ui.verbose
            shape_ui('PostTitle', toc(t0));
        end
        
        % -----------------------------------------------------------------
        % Update matching term
        if opt.ui.verbose
            shape_ui('Title', 'Update Warped Templates', false);
        end
        if ~isfield(opt.par.subjects.job, 'mem_warp')
            opt.par.subjects.job.mem_warp = opt.par.subjects.job.mem;
        end
        opt.par.subjects.job.mem = opt.par.subjects.job.mem_warp;
        [opt.par.subjects, dat] = distribute(opt.par.subjects, ...
            'updateWarpedTemplate', 'inplace', dat, model, opt);
        opt.par.subjects.job.mem_warp = opt.par.subjects.job.mem;
        if opt.ui.verbose
            shape_ui('Title', 'Update Matching terms', false);
        end
        if ~isfield(opt.par.subjects.job, 'mem_match')
            opt.par.subjects.job.mem_match = opt.par.subjects.job.mem;
        end
        opt.par.subjects.job.mem = opt.par.subjects.job.mem_match;
        [opt.par.subjects, dat] = distribute(opt.par.subjects, ...
            'lbMatching', 'inplace', dat, model, opt);
        opt.par.subjects.job.mem_match = opt.par.subjects.job.mem;
        model = aggregateMatching(dat, model, opt);
        model = updateLowerBound(model, opt);
        shape_plot_all(model,opt);
    end
    
    
    % =====================================================================
    % Write (if needed)
    dat   = structToFile(dat,   datpath);
    model = structToFile(model, modelpath);
    
    % =====================================================================
    % Save everything (to allow starting from a previous state)
    if ~isempty(opt.fnames.result)
        if opt.ui.verbose
            t0 = shape_ui('Title', 'Save model', false, true);
        end
        % Ensure nifti headers are ok
        createAllNifti(dat, model, opt);
        % Write workspace
        ftrack = opt.ui.ftrack;
        opt.ui.ftrack = nan;
        save(fullfile(opt.dir.model, opt.fnames.result), ...
             'model', 'dat', 'opt');
        opt.ui.ftrack = ftrack;
        if opt.ui.verbose
            shape_ui('PostTitle', toc(t0));
        end
    end
end