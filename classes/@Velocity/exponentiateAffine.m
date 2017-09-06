function A = exponentiateAffine(obj, q, B)
% FORMAT A = obj.exponentiateAffine((q), (B))
% q - Affine transform parameters in the Lie algebra
% B - Basis matrices of the ie algebra
% A - Reconstructed affine transform exp(sum_k q_k * B_k)

    % --- Check if nothing to do
    if nargin == 1 && obj.utd.A
        A = obj.A;
        return
    end

    % --- Default parameters
    if nargin < 3
        B = obj.AffineBasis;
        if nargin < 2
            q = obj.q;
        end
    end
    
    % --- Check input
    if ~obj.checkarray(B)
        if obj.Debug
            error('Cannot exponentiate affine transform. Missing basis')
        end
    end
    if obj.Debug,  fprintf('* exponentiateAffine\n'); end;
    
    % Compute Lie algebra matrix
    if isempty(q)
        A = eye(4);
    else
        q = reshape(q, [1 1 length(q)]);
        Q = sum(bsxfun(@times, B, q), 3);
        A = expm(Q);
    end
    
    % Store if needed
    if nargout == 0
        obj.A = A;
        obj.statusChanged('A');
    end