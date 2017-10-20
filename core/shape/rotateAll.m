function [model, dat] = rotateAll(model, dat, opt, R, iR)
% FORMAT [model, dat] = rotateAll(model, dat, opt, R, iR)

% "Rotate" the latent space.
%
% If R is orthogonal (true rotation), this causes no change at all to the
% model likelihood and the reconstructed velocities. The idea is then to
% correct for the lack of orthogonality in the PPCA framework (see Bishop's
% chapter). Because we make strong approximations when updating the
% principal subspace, this may improve the overall convergence.
% Additionally, it allows ro recover the true principal components (and not
% just a principal subspace).
% Here, iR is not needed.
%
% If R is not orthogonal, then iR must be provided. The idea is to combine
% a rotation, as described above, with a rescaling of the principal
% components (and corresponding inverse rescaling of the latent
% coordinates) in order to maximise p( z | W ) and p( W ).

    if opt.debug, fprintf(' * rotateAll\n'); end;
    
    if nargin < 5
        if norm(R*R' - eye(size(R, 1))) > 1e-7
            warning('R is not orthogonal')
        end
        iR = R';
    end
    
    % Rotate W
    for z=1:size(model.w,3)
        w1  = single(model.w(:,:,z,:,:));
        dim = [size(w1) 1 1 1];
        w1  = reshape(w1, [], opt.K) * iR;
        model.w(:,:,z,:,:) = reshape(w1, dim);
    end
    
    % Rotate sufficient statistics
    model.ww(:,:) = iR' * numeric(model.ww) * iR;
    model.z(:,:)  = R   * numeric(model.z);
    model.zz(:,:) = R   * numeric(model.zz) * R';
    model.Sz(:,:) = R   * numeric(model.Sz) * R';
    
    % Rotate subjects
    dat = batchProcess('Rotate', dat, opt, R);

end