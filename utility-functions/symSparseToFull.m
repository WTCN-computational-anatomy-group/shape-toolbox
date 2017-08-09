function f = symSparseToFull(s)
% FORMAT f = symSparseToFull(s)
%
% Convert a (field of) sparse symmetric squre matrix to a (field of) full
% squre matrix.
    
    % --- Dim info
    dim = size(s);
    lat = dim(1:end-1);
    k   = dim(end);
    [ind, n] = symIndices(k);
    
    % --- Allocate output
    s = reshape(numeric(s), [], k);
    f = zeros([prod(lat) n n], class(s));
    
    % --- Transfer data
    for l1=1:n
        for l2=1:n
            f(:,l1,l2) = s(:,ind(l1,l2));
        end
    end
    
    % --- Reshape to original
    f = reshape(f, [lat n n]);
end