function m = loadDiag(m)
% Additional regularisation in case a matrix is singular

    factor = 1e-7;
    while rcond(m) < 1e-5
        m = m + factor * max(diag(m)) * eye(size(m));
        factor = 10 * factor;
    end

end