function h = pointwiseInv(h, sym)
% Performs the pointwise inversion of a tensor field. The last two
% dimensions are assumed to be that of the matrices, and the other that of
% the lattice.

%     h = h(14,5,:,:);

    if nargin < 2
        sym = false;
    end

    % --- Dim info
    dim = size(h);
    if sym
        K = dim(end);
        lat = dim(1:end-1);
        [ind, N] = symIndices(K);
        M = N;
    else
        N = dim(end-1);
        M = dim(end);
        lat = dim(1:end-2);
        K = N*M;
    end
    
    % --- Helper function
    function X = invert(X)
        if sym
            Y = zeros(N,M);
            for n=1:N
                for m=1:M
                    Y(n,m) = X(ind(n,m));
                end
            end
            X = Y;
            clear Y;
        else
            X = reshape(X, N, M);
        end
        % In case X is not diagonal dominant
        blockdiag = false;
        if issame(X(:,end), [zeros(N-1,1) ; 1]) && issame(X(end,:), [zeros(1,M-1) 1])
            blockdiag = true;
            X = X(1:end-1,1:end-1);
        end
        R = max(0, max(sum(abs(X - diag(diag(X))), 2) - diag(abs(X))));
        if R > 0
            X = X + max(diag(X)) * eye(size(X));
        end
        X = inv(X);
        if blockdiag
            X = [X zeros(N-1,1); zeros(1,M-1) 1];
        end
    end

    % --- Pointwise operation
    h = numeric(h);
    h = num2cell(reshape(h, [], K), 2);
    h = cellfun(@invert, h, 'UniformOutput', false);
    h = cat(3, h{:});
    h = reshape(permute(h, [3 2 1]), [lat N M]);
end