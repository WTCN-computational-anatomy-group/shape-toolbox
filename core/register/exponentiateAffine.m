function A = exponentiateAffine(q, B, varargin)
% FORMAT A = exponentiateAffine(q, B, ...)
%
% ** Required **
% q - Affine transform parameters in the Lie algebra
% B - Basis matrices of the Lie algebra
% ** Output **
% A - Reconstructed affine transform = exp(sum_k q_k * B_k)

    % --- Parse inputs
    p = inputParser;
    p.FunctionName = 'exponentiateAffine';
    p.addRequired('q');
    p.addRequired('B', @checkarray);
    p.addParameter('output', []);
    p.addParameter('debug', false);
    p.parse(q, B, varargin{:});
    
    if p.Results.debug, fprintf('* exponentiateAffine\n'); end;
    
    if isempty(q)
        A = eye(4);
    else
        % --- Load data
        q = numeric(q);
        B = numeric(B);
        % --- Compute Lie algebra matrix
        q = reshape(q, [1 1 length(q)]);
        Q = sum(bsxfun(@times, B, q), 3);
        A = expm(Q);
    end
    
    % --- Write on disk
    if ~isempty(p.Results.output)
        A = saveOnDisk(p.Results.output, A, 'name', 'A');
    end
    
end