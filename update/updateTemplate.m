function model = updateTemplate(model, opt)
% FORMAT dat = updateSubspace(dat, model, opt)
% model - Model-specific data
% opt   - Options
%
% All inputs are structures that can either be in memory or on disk in the
% form of a mat file. In the latter case, it is read and, if needed,
% written back.
%--------------------------------------------------------------------------
% Update template by Gauss-Newton (smooth prior)
%__________________________________________________________________________
% Copyright (C) 2018 Wellcome Centre for Human Neuroimaging

    % =====================================================================
    % Read input from disk (if needed)
    [model, modelpath, opt] = fileToStruct(model, opt);
    
    % =====================================================================
    % Setup some constant
    spm_field('boundary', opt.tpl.bnd);
    
    % =====================================================================
    % /!\ In the categorical case, we have no choice but to load everything
    % in memory.
    if strcmpi(opt.model.name, 'categorical')
        
        % -----------------------------------------------------------------
        % Rotate out null space > force sum(a, 4) == 0
        a = single(numeric(model.tpl.a));
        a = rotateLogTemplate(...
            a, ...
            'par',   opt.par.within_main, ...
            'debug', opt.ui.debug);
        
        % -----------------------------------------------------------------
        % Add prior part of the gradient
        g = single(numeric(model.tpl.g));
        g = g + spm_field('vel2mom', single(a), double([opt.tpl.vs opt.tpl.prm]));
        h = single(numeric(model.tpl.h));
        
        % -----------------------------------------------------------------
        % Gauss-Newton update
        a = a - spm_field(h, g, double([opt.tpl.vs opt.tpl.prm 2 2]));
        clear h g
        
        % -----------------------------------------------------------------
        % Unrotate null space
        a = unrotateLogTemplate(...
            a, ...
            'par',   opt.par.within_main, ...
            'debug', opt.ui.debug);
        
        % -----------------------------------------------------------------
        % Write on disk
        model.tpl.a = saveOnDisk(model.tpl.a, a);
        
    % =====================================================================
    % In other cases, the Hessian is diagonal, so we can perform the
    % update one modality at a time.
    else
        
        for k=1:size(model.tpl.a,4)
            % -------------------------------------------------------------
            % Load one class/modality
            a = single(model.tpl.a(:,:,:,k));
            
            % -------------------------------------------------------------
            % Add prior part of the gradient
            g = single(model.tpl.g(:,:,:,k));
            g = g + spm_field('vel2mom', a, double([opt.tpl.vs opt.tpl.prm]));
            h = single(model.tpl.h(:,:,:,k));
            
            % -------------------------------------------------------------
            % Gauss-Newton update
            a = a - spm_field(h, g, double([opt.tpl.vs opt.tpl.prm 2 2]));
            
            % -------------------------------------------------------------
            % Write on disk
            model.tpl.a(:,:,:,k) = a;
        end
        
    end
    rmarray(model.tpl.g);
    rmarray(model.tpl.h);
    
    % =====================================================================
    % Exit
    model = structToFile(model, modelpath);
end