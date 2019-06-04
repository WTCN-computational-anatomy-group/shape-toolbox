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
% This function is (basically) pipelining all posterior and mode updates.
% It perform all updates the "right" (but less efficient) way to accurately 
% track the lower bound all along.
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
    %   AFFINE
    if model.q.active
        % -----------------------------------------------------------------
        % Update rigid/affine transform
        if opt.optimise.q.q
            if opt.ui.verbose
                shape_ui('Title', 'Update Affine', false);
            end
            okpre = sum(toArray(dat, '.q.ok') >= 0);
            if strcmpi(opt.par.subjects.mode, 'qsub')
                if ~isfield(opt.par.subjects.job, 'mem_affine')
                    opt.par.subjects.job.mem_affine = opt.par.subjects.job.mem;
                end
                opt.par.subjects.job.mem = opt.par.subjects.job.mem_affine;
            end
            [opt.par.subjects, dat] = distribute(opt.par.subjects, ...
                'updateAffine', 'inplace', dat, model, opt);
            if strcmpi(opt.par.subjects.mode, 'qsub')
                opt.par.subjects.job.mem_affine = opt.par.subjects.job.mem;
            end
            if opt.ui.verbose
                shape_ui('Title', '', false);
                okpost = sum(toArray(dat, '.q.ok') > 0);
                fprintf('success: %d / %d\n', okpost, okpre);
            end
            dat   = distribute('lbAffine', 'inplace', dat, model, opt);
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
            dat   = distribute('lbAffine', 'inplace', dat, model, opt);
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
            if strcmpi(opt.par.subjects.mode, 'qsub')
                if ~isfield(opt.par.subjects.job, 'mem_velocity')
                    opt.par.subjects.job.mem_velocity = opt.par.subjects.job.mem;
                end
                opt.par.subjects.job.mem = opt.par.subjects.job.mem_velocity;
            end
            [opt.par.subjects, dat] = distribute(opt.par.subjects, ...
                'updateVelocityShape', 'inplace', dat, model, opt);
            if strcmpi(opt.par.subjects.mode, 'qsub')
                opt.par.subjects.job.mem_velocity = opt.par.subjects.job.mem;
            end
            if opt.ui.verbose
                shape_ui('Title', '', false);
                okpost = sum(toArray(dat, '.v.ok') > 0);
                fprintf('success: %d / %d\n', okpost, okpre);
            end
            dat   = distribute('lbVelocityShape', 'inplace', dat, model, opt);
            model = aggregateVelocity(dat, model, opt);
            model = aggregateMatching(dat, model, opt);
            model = updateLowerBound(model, opt);
            shape_plot_all(model,opt);
        end
        
        % -----------------------------------------------------------------
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
            dat   = distribute('lbVelocityShape', 'inplace', dat, model, opt);
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
            dat   = distribute('updateResidualLB', 'inplace', dat, model, opt);
            dat   = distribute('lbLatent', 'inplace', dat, model, opt);
            dat   = distribute('lbVelocityShape', 'inplace', dat, model, opt);
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
            if strcmpi(opt.par.subjects.mode, 'qsub')
                if ~isfield(opt.par.subjects.job, 'mem_latent')
                    opt.par.subjects.job.mem_latent = opt.par.subjects.job.mem;
                end
                opt.par.subjects.job.mem = opt.par.subjects.job.mem_latent;
            end
            [opt.par.subjects, dat] = distribute(opt.par.subjects, ...
                'updateLatent', 'inplace', dat, model, opt);
            if strcmpi(opt.par.subjects.mode, 'qsub')
                opt.par.subjects.job.mem_latent = opt.par.subjects.job.mem;
            end
            dat   = distribute('updateResidualLB', 'inplace', dat, model, opt);
            dat   = distribute('lbLatent', 'inplace', dat, model, opt);
            dat   = distribute('lbVelocityShape', 'inplace', dat, model, opt);
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
            dat   = distribute('rotateLatent', Q, 'inplace', dat);
            dat   = distribute('lbLatent', 'inplace', dat, model, opt);
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
            dat   = distribute('lbLatent', 'inplace', dat, model, opt);
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
        if strcmpi(opt.par.subjects.mode, 'qsub')
            if ~isfield(opt.par.subjects.job, 'mem_template')
                opt.par.subjects.job.mem_template = opt.par.subjects.job.mem;
            end
            opt.par.subjects.job.mem = opt.par.subjects.job.mem_template;
        end
        [opt.par.subjects, dat] = distribute(opt.par.subjects, ...
            'updateTemplateGradHess', 'inplace', dat, model, opt);
        if strcmpi(opt.par.subjects.mode, 'qsub')
            opt.par.subjects.job.mem_template = opt.par.subjects.job.mem;
        end
            
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
        if strcmpi(opt.par.subjects.mode, 'qsub')
            if ~isfield(opt.par.subjects.job, 'mem_warp')
                opt.par.subjects.job.mem_warp = opt.par.subjects.job.mem;
            end
            opt.par.subjects.job.mem = opt.par.subjects.job.mem_warp;
        end
        [opt.par.subjects, dat] = distribute(opt.par.subjects, ...
            'updateWarpedTemplate', 'inplace', dat, model, opt);
        if strcmpi(opt.par.subjects.mode, 'qsub')
            opt.par.subjects.job.mem_warp = opt.par.subjects.job.mem;
        end
        if opt.ui.verbose
            shape_ui('Title', 'Update Matching terms', false);
        end
        if strcmpi(opt.par.subjects.mode, 'qsub')
            if ~isfield(opt.par.subjects.job, 'mem_match')
                opt.par.subjects.job.mem_match = opt.par.subjects.job.mem;
            end
            opt.par.subjects.job.mem = opt.par.subjects.job.mem_match;
        end
        [opt.par.subjects, dat] = distribute(opt.par.subjects, ...
            'lbMatching', 'inplace', dat, model, opt);
        if strcmpi(opt.par.subjects.mode, 'qsub')
            opt.par.subjects.job.mem_match = opt.par.subjects.job.mem;
        end
        model = aggregateMatching(dat, model, opt);
        model = updateLowerBound(model, opt);
        shape_plot_all(model,opt);
    end
    
    
    % =====================================================================
    % Write (if needed)
    dat   = structToFile(dat,   datpath);
    model = structToFile(model, modelpath);
end