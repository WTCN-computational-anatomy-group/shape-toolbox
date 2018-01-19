function [T, iT] = orthogonalisationMatrix(ezz, ww)
% FORMAT orthogonalisationMatrix(ezz, ww)
%
% ** Input **
% ezz - E[ZZ'] = E[Z]E[Z]' + cov[Z]
% ww  - Precision matrix for Z = W'*L*W
% ** Output **
% T   - Orthogonalisation matrix,  s.t.  T * ezz * T' ~ diag
% iT  - Almost inverse of T,       s.t. iT' * ww * iT ~ diag
%
% Copied from John's code. I still don't get everything.


    [Vz, Dz2]  = svd(double(ezz));
    [Vw, Dw2]  = svd(double(ww));
    Dz         = diag(sqrt(diag(Dz2) + eps));
    Dw         = diag(sqrt(diag(Dw2) + eps));
    [U, D, V]  = svd(Dw * Vw' * Vz * Dz');
    Dz         = spm_matcomp('LoadDiag', Dz);
    Dw         = spm_matcomp('LoadDiag', Dw);
    T          = D * V' * (Dz \ Vz');
    iT         = Vw * (Dw \ U);
    
end