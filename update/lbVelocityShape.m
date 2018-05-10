function dat = lbVelocityShape(dat, model, opt)
% FORMAT dat = lbVelocityShape(dat, model, opt)
% dat   - Subject-specific data
% model - Model-specific data
% opt   - Options
%
% All inputs are structures that can either be in memory or on disk in the
% form of a mat file. In the latter case, it is read and, if needed,
% written back.
%--------------------------------------------------------------------------
% Update velocity Log-likelihood/KL-divergence (dat.v.lb.val/ll1/ll2)
%__________________________________________________________________________
% Copyright (C) 2018 Wellcome Centre for Human Neuroimaging

    % =====================================================================
    % Read input from disk (if needed)
    [dat, datpath, model, ~, opt] = fileToStruct(dat, model, opt);

    % =====================================================================
    % Lower bound
    K = prod(opt.tpl.lat)*3;
    dat.v.lb.uncty = trace(dat.z.S*model.pg.ww);
    if defval(dat.v, '.observed', false)
        % -----------------------------------------------------------------
        % Observed velocity field
        
        % Shape part (Gaussian log-likelihood)
        dat.v.lb.ll1 = -0.5*( K*log(2*pi) ...
                              - K * model.v.LogLambda ...
                              - opt.pg.LogDetL ...
                              + model.v.l * dat.v.lb.reg ...
                              + model.v.l * dat.v.lb.uncty );
        % Shape-agnostic part (Gaussian log-likelihood)
        dat.v.lb.ll2 = -0.5*( K * log(2*pi) ...
                              - opt.pg.LogDetL ...
                              + dat.v.lb.regv );
        % Log-likelihood of a mixture of Gaussians
        dat.v.lb.val = model.mixreg.w(1) * dat.v.lb.ll1 + ...
                       model.mixreg.w(2) * dat.v.lb.ll2;
    else
        % -----------------------------------------------------------------
        % Latent velocity field
        
        % Shape part (Gaussian log-likelihood)
        dat.v.lb.ll1 = -0.5*( K*log(2*pi) ...
                              - K * model.v.LogLambda ...
                              - opt.pg.LogDetL ...
                              + model.v.l * dat.v.lb.reg ...
                              + model.v.l * dat.v.lb.uncty ...
                              + model.v.l * dat.v.lb.tr );
        % Shape-agnostic part (Gaussian log-likelihood)
        dat.v.lb.ll2 = -0.5*( K * log(2*pi) ...
                              - opt.pg.LogDetL ...
                              + dat.v.lb.regv ...
                              + dat.v.lb.tr );
        % KL-divergence between a Gaussian (Laplace approximation) and a
        % mixture of Gaussians
        dat.v.lb.val = -0.5*( - K ...
                              - model.mixreg.w(1) * K * model.v.LogLambda ...
                              - opt.pg.LogDetL  ...
                              + model.mixreg.w(1) * model.v.l * dat.v.lb.reg ...
                              + model.mixreg.w(2) * dat.v.lb.regv ...
                              + model.mixreg.w(1) * model.v.l * dat.v.lb.tr ...
                              + model.mixreg.w(2) * dat.v.lb.tr ...
                              + model.mixreg.w(1) * model.v.l * dat.v.lb.uncty );
                          
        % I do not include the LogDet of the posterior matrix for now
        % because the stochastic approximation has a huge variance and it
        % makes the lower bound plot completely useless.
        % I good potentially add it at the end (with more MC samples) to
        % approximate the model evidence.
        % > + dat.v.lb.ld
    end

    % =====================================================================
    % Exit
    dat = structToFile(dat, datpath);
        
end