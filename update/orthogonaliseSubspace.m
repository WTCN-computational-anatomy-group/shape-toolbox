function [model,Q] = orthogonaliseSubspace(model, opt)
% FORMAT [model,Q] = orthogonaliseSubspace(model, opt)
% model - Model-specific data
% opt   - Options
%
% All inputs are structures that can either be in memory or on disk in the
% form of a mat file. In the latter case, it is read and, if needed,
% written back.
%--------------------------------------------------------------------------
% Orthogonalise subspace
%__________________________________________________________________________
% Copyright (C) 2018 Wellcome Centre for Human Neuroimaging

    % =====================================================================
    % Read input from disk (if needed)
    [model, modelpath, opt] = fileToStruct(model, opt);
                 
    % ---------------------------------------------------------------------
    % Orthogonalise
    [U, iU] = orthogonalisationMatrix(model.z.zz+model.z.S, model.pg.ww);

    % ---------------------------------------------------------------------
    % Rescale
    [Q, iQ] = gnScalePG(iU' * model.pg.ww * iU * model.pg.n, ...
                        U   * model.z.zz  * U', ...
                        U   * model.z.S   * U', ...
                        opt.z.A0, opt.z.n0, model.pg.n);
    Q  = Q  * U;
    iQ = iU * iQ;

    % ---------------------------------------------------------------------
    % Rotate subspace
    for z=1:size(model.pg.w,3)
        w1  = single(model.pg.w(:,:,z,:,:));
        dim = [size(w1) 1 1 1];
        w1  = reshape(w1, [], opt.pg.K) * iQ;
        model.pg.w(:,:,z,:,:) = reshape(w1, dim);
    end
    
    % ---------------------------------------------------------------------
    % Rotate sufficient statistics
    model.pg.ww = iQ' * model.pg.ww * iQ;
    model.z.z   = Q   * model.z.z;
    model.z.zz  = Q   * model.z.zz  * Q';
    model.z.S   = Q   * model.z.S   * Q';
    model.z.Z   = Q   * model.z.Z;
    
    % ---------------------------------------------------------------------
    % Update latent prior
    model = updateLatentPrior(model, opt);
            
    % =====================================================================
    % Exit
    model = structToFile(model, modelpath);
end