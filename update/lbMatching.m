function dat = lbMatching(dat, ~, opt)
% FORMAT dat = lbLatent(dat, model, opt)
% dat   - Subject-specific data
% model - Model-specific data
% opt   - Options
%
% All inputs are structures that can either be in memory or on disk in the
% form of a mat file. In the latter case, it is read and, if needed,
% written back.
%--------------------------------------------------------------------------
% Update latent coordinates KL divergence (dat.z.lb.val)
%__________________________________________________________________________
% Copyright (C) 2018 Wellcome Centre for Human Neuroimaging

    % =================================================================
    % Read input from disk (if needed)
    [dat, datpath,  opt] = fileToStruct(dat, opt);

    % =================================================================
    % KL-divergence
    if isfield(dat, 'model')
        noisemodel = dat.model;
    else
        noisemodel = opt.model;
    end
    dat.f.lb.val = llMatching(...
        noisemodel, ...         % Matching term (Categorical/Normal/...)
        dat.tpl.wmu, ...        % Warped+Softmaxed+Rescaled template
        dat.f.f, ...            % Observed image (class responsibilities)
        'par', opt.par.within_subject, ... % parallelise processing? (usually no)
        'debug', opt.ui.debug); % Write debugging stuff? (usually no)

    % =================================================================
    % Exit
    dat = structToFile(dat, datpath);
        
end