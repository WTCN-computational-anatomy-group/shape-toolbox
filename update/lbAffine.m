function dat = lbAffine(dat, model, opt)
% FORMAT dat = lbAffine(dat, model, opt)
% dat   - Subject-specific data
% model - Model-specific data
% opt   - Options
%
% All inputs are structures that can either be in memory or on disk in the
% form of a mat file. In the latter case, it is read and, if needed,
% written back.
%--------------------------------------------------------------------------
% Update affine KL divergence (dat.q.lb.val)
%__________________________________________________________________________
% Copyright (C) 2018 Wellcome Centre for Human Neuroimaging


    if opt.q.Mr

        % =================================================================
        % Read input from disk (if needed)
        [dat, datpath, model, ~, opt] = fileToStruct(dat, model, opt);


        % =================================================================
        % KL-divergence
        rind = opt.q.rind;
        qq   = dat.q.qq(rind,rind);
        Sq   = dat.q.S(rind,rind);
        dat.q.lb.type = 'kl';
        dat.q.lb.val  = -0.5*( trace((Sq + qq) * model.q.A) ...
                               - model.q.LogDetA ...
                               + spm_matcomp('LogDet', Sq) ...
                               - opt.q.Mr );
        dat.q.lb.type = 'kl';
        dat.q.lb.name = '-KL Affine';
    
        % =================================================================
        % Exit
    dat = structToFile(dat, datpath);
        
    end
end