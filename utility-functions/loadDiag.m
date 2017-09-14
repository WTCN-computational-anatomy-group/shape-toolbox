function m = loadDiag(m)
% Additional regularisation in case a matrix is singular

    while rcond(m) < 1e-5
        m = m + 1e-7 * max(diag(m)) * eye(size(m));
    end

end