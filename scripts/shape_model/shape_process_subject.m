function dat = shape_process_subject(dat, model, opt)
% FORMAT dat = shape_process_subject(dat, model, opt)
% dat   - Subject-specific data
% model - Model-specific data
% opt   - Options
%
% All inputs are structures that can either be in memory or on disk in the
% form of a mat file. In the latter case, it is read and, if needed,
% written back.
%--------------------------------------------------------------------------
% Perform subject-specific updates.
% This function can be run as an independent job on a cluster.
%__________________________________________________________________________
% Copyright (C) 2018 Wellcome Centre for Human Neuroimaging

    % =====================================================================
    % Read input from disk (if needed)
    [dat, datpath, model, ~, opt] = fileToStruct(dat, model, opt);
    
    if model.converged
        if opt.verbose
            fprintf('Already converged. Unnecessary iteration.\n');
        end
        return
    end
    
    % =====================================================================
    % Update lower bound due to "population" updates
    if opt.optimise.q.A
        dat = lbAffine(dat, model, opt);
    end
    if opt.optimise.z.A || opt.optimise.pg.w
        dat = lbLatent(dat, model, opt);
    end
    if opt.optimise.v.l
        dat = lbVelocityShape(dat, model, opt);
    end
    if opt.optimise.tpl.a
        dat = updateWarpedTemplate(dat, model, opt); % Pull and reconstruct
        dat = lbMatching(dat, model, opt);           % Update matching term
    end
    
    % =====================================================================
    % Update rigid/affine transform
    if opt.optimise.q.q
        dat = updateAffine(dat, model, opt);
        dat = lbAffine(dat, model, opt);
    end
    
    % =====================================================================
    % Update velocity
    if opt.optimise.v.v
        dat = updateVelocityShape(dat, model, opt);
        dat = lbVelocityShape(dat, model, opt);
    end
    
    % =====================================================================
    % Update latent coordinates
    if opt.optimise.z.z && model.emit > 1
        dat = updateLatent(dat, model, opt);
        dat = lbLatent(dat, model, opt);
        dat = lbVelocityShape(dat, model, opt); % < depends on E[zz']
    end
    
    % =====================================================================
    % [TODO] Update TPM scaling factors
%     if opt.optimise.tpl.scale
%         dat = updateTissueScale(dat, model, opt);
%     end
    
    % =====================================================================
    % Update template gradient/Hessian
    if opt.optimise.tpl.a && strcmpi(opt.tpl.update, 'map')
        dat = updateTemplateGradHess(dat, model, opt);
    end
    
    % =====================================================================
    % Exit
    dat = structToFile(dat, datpath);
end
