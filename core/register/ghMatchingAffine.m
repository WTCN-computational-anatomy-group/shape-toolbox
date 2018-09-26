function [g, h] = ghMatchingAffine(model, mu, f, gmu, A, B, varargin)
%__________________________________________________________________________
%
% Compute gradient/hessian of the **negative** log-likelihood with respect 
% to affine parameters.
%
% -------------------------------------------------------------------------
%
% FORMAT [(g),(h)] = ghMatchingAffine(model, mu, f, ga, A, B, ...)
% 
% REQUIRED
% --------
% model - Structure with fields:
%           * 'name'    : 'normal', 'laplace', 'bernoulli' or 'categorical'
%           * ('sigma2'): Normal variance  [1]
%           * ('b')     : Laplace variance [1]
% mu    - Reconstructed probability template [mx my mz nc]
% f     - Observed image [mx my mz nc]
% ga    - Template coefficients spatial gradients
%         > If ipsi provided, [mx my mz nc 3]
%         > Else              [nx ny nz nc 3]
% A     - Affine transform
% B     - Basis system of the affine Lie algebra
%
% OPTIONAL
% --------
% phi   - Direct (image-to-template) diffeomorphism [default: none]
% jac   - Jacobian of phi [default: recompute]
% 
% KEYWORD ARGUMENTS
% -----------------
% ipsi    - Inverse (subj to template) warp [mx my mz 3]
%           > If provided, compute "pushed(gradient(pulled))"
% circ    - (Push only) Boundary conditions for pushing [1]
% count   - Pushed voxel count ([nx ny nz nc])
%           > If provided, compute "gradient(pushed)"
% bb      - Bounding box (if different between template and pushed image)
% hessian - Only compute hessian (not gradient)
% Mmu     - Template voxel-to-world affine transform [default: eye(4)]
% approx  - Approximate hessian [default: true]
% loop    - How to split: 'none', 'component', 'slice' [default: auto]
% par     - Parallelise: false/true/number of workers [default: false]
%
% OUTPUT
% ------
% g     - First derivatives w.r.t. affine parameters ([nq])
% h     - Second derivatives w.r.t. affine parameters ([nq nq])
%__________________________________________________________________________

    % --- Parse inputs
    p = inputParser;
    p.FunctionName = 'ghMatchingAffine';
    p.addRequired('model',  @(X) isstruct(X) && isfield(X, 'name'));
    p.addRequired('mu',     @checkarray);
    p.addRequired('f',      @checkarray);
    p.addRequired('gmu',    @checkarray);
    p.addRequired('A',      @(X) issame(size(X), [4 4]));
    p.addRequired('basis',  @(X) issame(size(X(:,:,1)), [4 4]));
    p.addOptional('phi',       []);
    p.addOptional('jac',       []);
    p.addParameter('ipsi',     [],      @checkarray);
    p.addParameter('circ',     1);
    p.addParameter('count',    [],      @checkarray);
    p.addParameter('Mmu',      eye(4),  @(X) issame(size(X), [4 4]));
    p.addParameter('approx',   true,    @islogical);
    p.addParameter('bb',       struct,  @isstruct);
    p.addParameter('hessian',  false,   @islogical);
    p.addParameter('loop',     '',      @ischar);
    p.addParameter('par',      false,   @isscalar);
    p.addParameter('debug',    false,   @isscalar);
    p.parse(model, mu, f, gmu, A, B, varargin{:});
    phi     = p.Results.phi;
    jac     = p.Results.jac;
    c       = p.Results.count;
    ipsi    = p.Results.ipsi;
    circ    = p.Results.circ;
    Mmu     = p.Results.Mmu;
    bb      = p.Results.bb;
    approx  = p.Results.approx;
    hessian = p.Results.hessian;
    loop    = p.Results.loop;
    par     = p.Results.par;
    debug   = p.Results.debug;
    
    if debug, fprintf('* ghMatchingAffine\n'); end
    
    % --- Default bounding box
    dim = [size(mu) 1];
    if ~isfield(bb, 'x')
        bb.x = 1:dim(1);
    end
    if ~isfield(bb, 'y')
        bb.y = 1:dim(2);
    end
    if ~isfield(bb, 'z')
        bb.z = 1:dim(3);
    end
    bbx = bb.x;
    bby = bb.y;
    bbz = bb.z;

    % --- Gradient/Hessian w.r.t. initial velocity
    if hessian
        hv = ghMatchingVel(model, mu, f, gmu, ...
            'ipsi',     ipsi, ...
            'circ',     circ, ...
            'count',    c, ...
            'bb',       bb, ...
            'hessian',  true, ...
            'loop',     loop, ...
            'par',      par, ...
            'debug',    debug);
        hv = hv(bbx,bby,bbz,:);
    elseif nargin > 1
        [gv, hv] = ghMatchingVel(model, mu, f, gmu, ...
            'ipsi',     ipsi, ...
            'circ',     circ, ...
            'count',    c, ...
            'bb',       bb, ...
            'loop',     loop, ...
            'par',      par, ...
            'debug',    debug);
        gv = gv(bbx,bby,bbz,:);
        hv = hv(bbx,bby,bbz,:);
    else
        gv = ghMatchingVel(model, mu, f, gmu, ...
            'ipsi',     ipsi, ...
            'circ',     circ, ...
            'count',    c, ...
            'bb',       bb, ...
            'loop',     loop, ...
            'par',      par, ...
            'debug',    debug);
        gv = gv(bbx,bby,bbz,:);
    end
    
    % --- Check diffeomorphic part
    jac = single(numeric(jac));
    phi = single(numeric(phi));
    if isempty(phi)
        id = spm_warps('identity', dim(1:3));
    elseif isempty(jac)
        jac = spm_diffeo('def2jac', phi);
    end
    jac = jac(bbx,bby,bbz,:,:);
    
    % --- Allocate gradient and hessian
    do_gradient = ~hessian;
    do_hessian  = nargout > 1 || hessian;
    nq = size(B, 3);
    if do_gradient
        g = zeros([nq 1]);
    end
    if do_hessian
        h = zeros(nq);
    end
    
    % --- Gradient + Symmetric Hessian part
    for i=1:nq
        dXi = Mmu \ A \ B(:,:,i) * A * Mmu;
        if ~isempty(phi)
            dXi = spm_warps('transform', dXi, phi);
            dXi = dXi(bbx,bby,bbz,:);
            dXi = spm_matcomp('Pointwise3', jac, dXi, 'i');
        else
            dXi = spm_warps('transform', dXi, id);
            dXi = dXi(bbx,bby,bbz,:);
        end
        if do_gradient
            g(i) = sumall(spm_matcomp('Pointwise3', gv, dXi));
        end
        if do_hessian
            for j=1:i
                dXj = Mmu \ A \ B(:,:,j) * A * Mmu;
                if ~isempty(phi)
                    dXj = spm_warps('compose', dXj, phi);
                    dXj = dXj(bbx,bby,bbz,:);
                    dXj = spm_matcomp('Pointwise3', jac, dXj, 'i');
                else
                    dXj = spm_warps('compose', dXj, id);
                    dXj = dXj(bbx,bby,bbz,:);
                end
                h(i,j) = h(i,j) + sumall(spm_matcomp('Pointwise3', dXi, spm_matcomp('Pointwise3', hv, dXj)));
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
    if do_hessian && ~approx
        for i=1:nq
            for j=1:nq
                Bij = B(:,:,i) * B(:,:,j);
                if any(any(Bij))
                    dXij =  Mmu \ A \ Bij * A * Mmu;
                    if ~isempty(phi)
                        dXij = spm_warps('compose', dXij, phi);
                        dXij = dXij(bbx,bby,bbz,:);
                        dXij = spm_matcomp('Pointwise3', jac, dXij, 'i');
                    else
                        dXij = spm_warps('compose', dXij, id);
                        dXij = dXij(bbx,bby,bbz,:);
                    end
                    h(i,j) = h(i,j) + sumall(spm_matcomp('Pointwise3', gv, dXij));
                end
            end
        end
        h = (h + h')/2;   % Insure symmetric Hessian
    end

    if hessian
        g = [];
        [g, h] = deal(h, g);
    end
    
end
    