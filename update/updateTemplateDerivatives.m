function model = updateTemplateDerivatives(model, opt)
% FORMAT model = updateTemplateDerivatives(model, opt)
% model - Model-specific data
% opt   - Options
%
% All inputs are structures that can either be in memory or on disk in the
% form of a mat file. In the latter case, it is read and, if needed,
% written back.
%--------------------------------------------------------------------------
% Update spatial gradients + reconstruct template (for visualisation)
%__________________________________________________________________________
% Copyright (C) 2018 Wellcome Centre for Human Neuroimaging

    % =====================================================================
    % Read input from disk (if needed)
    [model, modelpath, opt] = fileToStruct(model, opt);
                 
    % =====================================================================
    % Compute spatial derivatives
    model.tpl.gmu = templateGrad(...
        model.tpl.a,  ...           % (Log)-template coefficients
        ~opt.tpl.bnd,  ...          % Boundary conditions
        'debug',  opt.ui.debug, ... % Write debugging stuff?
        'output', model.tpl.gmu);   % Output file array (saves memory)
    
    % =====================================================================
    % 2D case: ensure null derivatives in 3rd direction
    if opt.model.dim == 2
        for z=1:size(model.tpl.gmu, 3)
            model.tpl.gmu(:,:,z,:,3) = 0;
        end
    end
    
    % =====================================================================
    % Reconstruct probability template
    if opt.tpl.cat
        model.tpl.mu = reconstructProbaTemplate(...
            model.tpl.a, ...               % (Log)-template coefficients
            'par',    opt.par.within_main, ... % Parallelise processing?
            'debug',  opt.ui.debug,  ...   % Write debugging stuff?
            'output', model.tpl.mu);       % Output file array (saves memory)
    else
        model.tpl.mu = model.tpl.a;
    end

    % =====================================================================
    % Exit
    model = structToFile(model, modelpath);
end