function dat = lbVelocity(dat, ~, opt)
% FORMAT dat = lbVelocity(dat, model, opt)
% dat   - Subject-specific data
% model - Model-specific data
% opt   - Options
%
% All inputs are structures that can either be in memory or on disk in the
% form of a mat file. In the latter case, it is read and, if needed,
% written back.
%--------------------------------------------------------------------------
% Update velocity Log-likelihood/KL-divergence (dat.v.lb.val)
%__________________________________________________________________________
% Copyright (C) 2018 Wellcome Centre for Human Neuroimaging

    % =====================================================================
    % Read input from disk (if needed)
    [dat, datpath, opt] = fileToStruct(dat, opt);

    % =====================================================================
    % Lower bound
    K = prod(opt.tpl.lat)*3;
    switch lower(opt.v.update)
        
        case 'mode'
            
            dat.v.lb.val = -0.5*( K * log(2*pi) ...
                                  - opt.v.LogDetL ...
                                  + dat.v.lb.reg );
            dat.v.lb.type = 'll';
            dat.v.lb.name = 'Velocity prior';
            
        case 'variational'
            
            dat.v.lb.val = -0.5*( - K ...
                                  - opt.v.LogDetL  ...
                                  + dat.v.lb.ld ...
                                  + dat.v.lb.reg ...
                                  + dat.v.lb.tr );
            dat.v.lb.type = 'kl';
            dat.v.lb.name = '-KL Velocity';
            
    end

    % =====================================================================
    % Exit
    dat = structToFile(dat, datpath);
        
end