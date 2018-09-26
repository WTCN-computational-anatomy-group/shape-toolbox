function dat = initAffineLaplace(dat, model, opt)
% FORMAT dat = initAffineLaplace(dat, model, opt)
% dat   - Subject-specific data
% model - Model-specific data
% opt   - Options
%
% All inputs are structures that can either be in memory or on disk in the
% form of a mat file. In the latter case, it is read and, if needed,
% written back.
%--------------------------------------------------------------------------
% Initialise affine Laplace approximation (i.e., compute Hessian)
%__________________________________________________________________________
% Copyright (C) 2018 Wellcome Centre for Human Neuroimaging

    % =====================================================================
    % Read input from disk (if needed)
    [dat, datpath, model, ~, opt] = fileToStruct(dat, model, opt);

    % =====================================================================
    % If the velocity is an observed -> nothing to do
    if defval(dat.v, '.observed', false)
        dat = structToFile(dat, datpath);
        return
    end
    
    % =====================================================================
    % Update Hessian for Laplace approximation
    if opt.q.Mr

        % Likelihood part
        h = ghMatchingAffine(...
            noisemodel, ...             % Matching model (categorical/normal/...)
            dat.tpl.wmu, ...            % Warped (+ softmaxed) template
            dat.f.f, ...                % Observed matched image (responsibility)
            model.tpl.gmu, ...          % (Log)-template spatial gradients
            dat.q.A, ...                % Current rigid/affine transform
            opt.q.B, ...                % Rigid/affine Lie basis
            phi, ...                    % Direct diffeomorphism
            jac, ...                    % Jacobian determinant of the direct diffeo
            'ipsi',    dat.v.ipsi, ...  % Complete (rigid+diffeo) inverse transform
            'hessian', true, ...        % Do not compute gradient
            'Mmu',     model.tpl.M, ... % Template voxel-to-world
            'approx',  opt.q.hapx, .... % Approximate hessian? (usually true)
            'par',     opt.par.within_subject, ... % Parallelise stuff? (usually no)
            'debug',   opt.ui.debug);   % Write debuging stuff? (usually no)

        % prior part
        [~, hq] = ghPriorAffine(dat.q(rind), model.q.A, 'debug', opt.ui.debug);
        h(rind,rind) = h(rind,rind) + hq;
        clear hq

        % Additional regularisation for robustness
        h = spm_matcomp('LoadDiag', h);
        
        % Compute covariance
        dat.q.S = spm_matcomp('Inv', h);
        
    else
        dat.q.S = 0;
    end
              
    % =====================================================================
    % Exit
    dat = structToFile(dat, datpath);
end