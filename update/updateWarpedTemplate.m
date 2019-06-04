function dat = updateWarpedTemplate(dat, model, opt)
% FORMAT dat = updateWarpedTemplate(dat, model, opt)
% dat   - Subject-specific data
% model - Model-specific data
% opt   - Options
%
% All inputs are structures that can either be in memory or on disk in the
% form of a mat file. In the latter case, it is read and, if needed,
% written back.
%--------------------------------------------------------------------------
% Update warped template (if template/deformation/scales changed)
%__________________________________________________________________________
% Copyright (C) 2018 Wellcome Centre for Human Neuroimaging

    % =====================================================================
    % Read input from disk (if needed)
    [dat, datpath, model, ~, opt, ~] = fileToStruct(dat, model, opt);
    
    % =====================================================================
    % If the velocity is an observed -> nothing to do
    if defval(dat.v, '.observed', false)
        dat = structToFile(dat, datpath);
        return
    end
                 
    % =====================================================================
    % Pull (log)-template
    if opt.optimise.tpl.a
        if opt.tpl.cat, wout = dat.tpl.wa;
        else,           wout = dat.tpl.wmu; end
        wout = pullTemplate(...
            dat.v.ipsi, ...             % Inverse diffeo+rigid transform
            model.tpl.a, ...            % Template (log)-coefficients
            'order',  opt.tpl.itrp, ... % Interpolation order
            'bnd',   ~opt.tpl.bnd, ...  % Boundary condition
            'output', wout, ...         % Output file_array (saves memory)
            'par',    opt.par.within_subject, ... % parallelise processing? (usually no)
            'debug',  opt.ui.debug);    % Write debugging stuff? (usually no)
        if opt.tpl.cat, dat.tpl.wa  = wout;
        else,           dat.tpl.wmu = wout; end
    end
    
    % =====================================================================
    % Softmax + rescale (if categorical)
    if opt.tpl.cat && (opt.optimise.tpl.a || opt.optimise.tpl.scale)
        dat.tpl.wmu = reconstructProbaTemplate(...
            dat.tpl.wa, ...               % Warped log-template
            'scale',  dat.tpl.scale, ...  % Rescaling factors
            'output', dat.tpl.wmu, ...    % Output file_array
            'par',    opt.par.within_subject, ... % parallelise processing? (usually no)
            'debug',  opt.ui.debug);      % Write debugging stuff? (usually no)
        if opt.optimise.tpl.a
            dat.tpl.wa = rmarray(dat.tpl.wa); % Remove warped-coefficients (saves disk)
        end
    end

    % =====================================================================
    % Exit
    dat = structToFile(dat, datpath);
end