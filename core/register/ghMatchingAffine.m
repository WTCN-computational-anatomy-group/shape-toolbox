function [g, h] = ghMatchingAffine(model, mu, f, c, gmu, A, B, varargin)
% FORMAT [g, (h)] = ghMatchingAffine(mu, f, c, gmu, A, basis, (phi), (jac), ...)
% 
% ** Required **
% model - 'normal', 'laplace', 'bernoulli', 'categorical'
% mu    - (Reconstructed) template image ([nx ny nz nc])
% f     - Observed image pushed into template space ([nx ny nz nc])
% c     - Pushed voxel count ([nx ny nz nc])
% gmu   - Template spatial gradients.
% A     - Affine transform
% B     - Basis system of the affine Lie algebra
% ** Optional **
% phi   - Direct (image-to-template) diffeomorphism [default: none]
% jac   - Jacobian of phi [default: recompute]
% ** Keyword arguments **
% sigma2- Noise variance for the normal case [default: 1]
% b     - Noise variance for the laplace case [default: 1]
% Mmu   - Template voxel-to-world affine transform [default: eye(4)]
% approx- Approximate hessian [default: true]
% loop  - How to split: 'none', 'component', 'slice' [default: auto]
% par   - Parallelise: false/true/number of workers [default: false]
% ** Outputs **
% g     - First derivatives w.r.t. affine parameters ([nq])
% h     - Second derivatives w.r.t. affine parameters ([nq nq])
%
% Compute gradient/hessian with respect to affine parameters.

    % --- Parse inputs
    p = inputParser;
    p.FunctionName = 'ghMatchingAffine';
    p.addRequired('model',  @ischar);
    p.addRequired('mu',     @checkarray);
    p.addRequired('f',      @checkarray);
    p.addRequired('c',      @checkarray);
    p.addRequired('gmu',    @checkarray);
    p.addRequired('A',      @(X) issame(size(X), [4 4]));
    p.addRequired('basis',  @(X) issame(size(X(:,:,1)), [4 4]));
    p.addOptional('phi',       []);
    p.addOptional('jac',       []);
    p.addParameter('sigma2',   1,       @isscalar);
    p.addParameter('b',        1,       @isscalar);
    p.addParameter('Mmu',      eye(4),  @(X) issame(size(X), [4 4]));
    p.addParameter('approx',   true,    @islogical);
    p.addParameter('loop',     '',      @ischar);
    p.addParameter('par',      false,   @isscalar);
    p.addParameter('output',   []);
    p.addParameter('debug',    false,   @isscalar);
    p.parse(model, mu, f, c, gmu, A, B, varargin{:});
    phi = p.Results.phi;
    jac = p.Results.jac;
    Mmu = p.Results.Mmu;
    
    if p.Results.debug, fprintf('* ghMatchingAffine\n'); end;
    
    % --- Gradient/Hessian w.r.t. initial velocity
    [gv, hv] = ghMatchingVel(model, mu, f, c, gmu, ...
        'sigma2',   p.Results.sigma2, ...
        'b',        p.Results.b, ...
        'loop',     p.Results.loop, ...
        'par',      p.Results.par, ...
        'debug',    p.Results.debug);
    
    % --- Check diffeomorphic part
    jac = single(numeric(jac));
    phi = single(numeric(phi));
    if isempty(phi)
        lat = [size(mu) 1];
        lat = lat(1:3);
        id = warps('identity', lat);
    elseif isempty(jac)
        jac = spm_diffeo('def2jac', phi);
    end
    
    % --- Allocate gradient and hessian
    do_hessian = nargout > 1;
    nq = size(B, 3);
    g = zeros([nq 1]);
    if do_hessian
        h = zeros(nq);
    end
    
    % --- Gradient + Symmetric Hessian part
    for i=1:nq
        dXi = Mmu \ A \ B(:,:,i) * A * Mmu;
        if ~isempty(phi)
            dXi = warps('transform', dXi, phi);
            dXi = pointwise3(jac, dXi, 'i');
        else
            dXi = warps('transform', dXi, id);
        end
        g(i) = sumall(pointwise3(gv, dXi));
        if do_hessian
            for j=1:i
                dXj = Mmu \ A \ B(:,:,j) * A * Mmu;
                if ~isempty(phi)
                    dXj = warps('compose', dXj, phi);
                    dXj = pointwise3(jac, dXj, 'i');
                else
                    dXj = warps('compose', dXj, id);
                end
                h(i,j) = h(i,j) + sumall(pointwise3(dXi, pointwise3(hv, dXj)));
            end
        end
    end
    if do_hessian
        for i=1:nq
            for j=(i+1):nq
                h(i,j) = h(j,i);
            end
        end
    end
    
    % --- Nonsymmetric Hessian part
    if do_hessian && ~p.Results.approx
        for i=1:nq
            for j=1:nq
                Bij = B(:,:,i) * B(:,:,j);
                if any(any(Bij))
                    dXij =  obj.Mmu \ A \ Bij * A * obj.Mmu;
                    if ~isempty(phi)
                        dXij = warps('compose', dXij, phi);
                        dXij = pointwise3(jac, dXij, 'i');
                    else
                        dXij = warps('compose', dXij, id);
                    end
                    h(i,j) = h(i,j) + sumall(pointwise3(gv, dXij));
                end
            end
        end
        h = (h + h')/2;   % Insure symmetric Hessian
    end
    
end
    