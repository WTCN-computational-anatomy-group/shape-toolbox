function [dat,model] = shape_process_pop(dat, model, opt)
% FORMAT [dat,model] = shape_process_pop(dat, model, opt)
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
% Nothing is distributed.
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
        end
    end
    
    % =====================================================================
    %   VELOCITY
    if model.v.active
        
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
            [~, dat] = distribute([], 'lbVelocityShape', 'inplace', dat, model, opt);
            model = aggregateVelocity(dat, model, opt);
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
            model = aggregateLatent(dat,model,opt);
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
        end
    end
    
    % =====================================================================
    %   TEMPLATE
    if opt.f.N && opt.optimise.tpl.a
            
        % -----------------------------------------------------------------
        % Update template
        if opt.ui.verbose
            t0 = shape_ui('Title', 'Update Template', false, true);
        end
        model = aggregateMatching(dat, model, opt);
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
    end
    
    model = updateLowerBound(model, opt);
    shape_plot_all(model,opt);
    
    % =====================================================================
    % Write (if needed)
    dat   = structToFile(dat,   datpath);
    model = structToFile(model, modelpath);
end