function dat = updateLatent(dat, model, opt)
% FORMAT dat = updateAffine(dat, model, opt)
% dat   - Subject-specific data
% model - Model-specific data
% opt   - Options
%
% All inputs are structures that can either be in memory or on disk in the
% form of a mat file. In the latter case, it is read and, if needed,
% written back.
%--------------------------------------------------------------------------
% Update latent shape coordinates
%__________________________________________________________________________
% Copyright (C) 2018 Wellcome Centre for Human Neuroimaging

    % =====================================================================
    % Read input from disk (if needed)
    [dat, datpath, model, ~, opt] = fileToStruct(dat, model, opt);
    
    % The posterior is a multivariate Gaussian with
    % covariance: S = w * lam * W'LW + A
    % mean:       z = w * lam * S * W'Lv
    
    % =====================================================================
    % Load stuff if needed
    if isempty(defval(dat.buffer, 'v.v', []))
        if opt.buf, dat.buffer.v.v = numeric(dat.v.v);
        else,       dat.buffer.v.v = dat.v.v;                   end
    end
    
    % =====================================================================
    % Update covariance
    dat.z.S = model.mixreg.w(1) * model.v.l * model.pg.ww + model.z.A;
    dat.z.S = spm_matcomp('Inv', dat.z.S);
   
    % =====================================================================
    % Update mean
    spm_diffeo('boundary', opt.pg.bnd);
    wm = zeros([opt.pg.K 1]);
    m = spm_diffeo('vel2mom', single(numeric(dat.buffer.v.v)), double([opt.tpl.vs opt.pg.prm]));
    for k=1:opt.pg.K
        w1    = model.pg.w(:,:,:,:,k);
        wm(k) = w1(:)'*m(:);
    end
    clear m w1
    dat.z.z  = model.mixreg.w(1) * model.v.l * dat.z.S * wm;
    dat.z.zz = dat.z.z * dat.z.z';
                             
    % =====================================================================
    % Exit
    dat = structToFile(dat, datpath);
end