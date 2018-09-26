function model = lbTemplate(model, opt)
% FORMAT model = lbTemplate(model, opt)
% model - Model-specific data
% opt   - Options
%
% All inputs are structures that can either be in memory or on disk in the
% form of a mat file. In the latter case, it is read and, if needed,
% written back.
%--------------------------------------------------------------------------
% Update prior log-likelihood for the (log)-template
%__________________________________________________________________________
% Copyright (C) 2018 Wellcome Centre for Human Neuroimaging

    % =====================================================================
    % Read input from disk (if needed)
    [model, modelpath, opt] = fileToStruct(model, opt);
     
    % =====================================================================
    % Lower Bound
    model.lb.a.val = 0;
    K = size(model.tpl.a,4);
    for k=1:K
        a1 = single(model.tpl.a(:,:,:,k));
        m1 = spm_field('vel2mom', a1, double([opt.tpl.vs opt.tpl.prm]));
        model.lb.a.val = model.lb.a.val - 0.5 * a1(:)'*m1(:);
        clear a1 m1
    end
    model.lb.a.val = model.lb.a.val ...
        - 0.5 * K*prod(opt.tpl.vs)*log(2*pi) ...
        + 0.5 * K*opt.tpl.LogDetL;
    model.lb.a.type = 'll';
    model.lb.a.name = 'Template prior';
                   
    % =====================================================================
    % Exit
    model = structToFile(model, modelpath);
end