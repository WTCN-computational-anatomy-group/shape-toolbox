function [g, h] = gradHessMatchingAffine(obj, mu, f, c, gmu, A, phi, jac)
% FORMAT [g, h] = obj.gradHessMatchingAffine((mu), (f), (c), (gmu))
% obj   - Velocity object
% (mu)  - Template image ([nx ny nz nc]) [default: obj.mu]
% (f)   - Observed image pushed into template space ([nx ny nz nc])
%         [default: obj.pf]
% (c)   - Pushed voxel count ([nx ny nz nc]) [default: obj.pvox]
% (gmu) - Template spatial gradients. [default: obj.gmu]
% g     - First derivatives w.r.t. full velocity ([nx ny nz nc])
% h     - Second derivatives w.r.t. full velocity
%         (diagonal approximation: [nx ny nz nc], except for the multinomial 
%         case where it is a symmetric tensor field of virtual size 
%         [nx ny nz nc nc])
%
% obj.MatchingTerm:
% normal      - Gaussian noise model     (F = Mu(IPhi) + N(0,s))
% laplace     - Laplace noise model      (F = Mu(IPhi) + L(0,b))
% binomial    - True/False realisation   (F = Ber(Mu(IPhi)))
% multinomial - Multiclass realisation   (F = Cat(Mu(IPhi)))
%
% Compute gradient/hessian with respect to affine parameters


    % --- Default arguments
    if nargin < 7
        obj.exponentiateVelocity('phi', 'jac');
        phi = obj.phi;
        jac = obj.jac;
    elseif nargin < 8
        obj.exponentiateVelocity('jac');
        jac = obj.jac;
    end
    if nargin < 6
        obj.exponentiateAffine();
        A = obj.A;
        if nargin < 5
            obj.computeTemplateGrad();
            gmu = obj.gmu;
            if nargin == 1
                obj.reconstructTemplate();
                mu = obj.mu;
                obj.computeTemplateGrad();
                gmu = obj.gmu;
                obj.pushImage();
                c = obj.pvox;
                f = obj.pf;
            elseif nargin < 4
                obj.pushImage();
                c = obj.pvox;
                if nargin < 3
                    f = obj.pf;
                end
            end
        end
    end
    
    [gv, hv] = obj.gradHessMatchingVel(mu, f, c, gmu);
    
    jac = single(numeric(jac));
    phi = single(numeric(phi));
    B   = obj.AffineBasis;
    
    if isempty(phi)
        lat = [size(mu) 1];
        lat = lat(1:3);
        id = warps('identity', lat);
    end
    
    nq = size(obj.AffineBasis, 3);
    g = zeros([nq 1]);
    h = zeros(nq);
    
    % --- Gradient + Symmetric Hessian part
    for i=1:nq
        dXi = obj.Mmu \ A \ B(:,:,i) * A * obj.Mmu;
        if ~isempty(phi)
            dXi = warps('transform', dXi, phi);
            dXi = pointwise3(jac, dXi, 'i');
        else
            dXi = warps('transform', dXi, id);
        end
        g(i) = sumall(pointwise3(gv, dXi));
        for j=1:i
            dXj = obj.Mmu \ A \ B(:,:,j) * A * obj.Mmu;
            if ~isempty(phi)
                dXj = warps('compose', dXj, phi);
                dXj = pointwise3(jac, dXj, 'i');
            else
                dXj = warps('compose', dXj, id);
            end
            h(i,j) = h(i,j) + sumall(pointwise3(dXi, pointwise3(hv, dXj)));
            clear dXj
        end
        clear dXi
    end
    for i=1:nq
        for j=(i+1):nq
            h(i,j) = h(j,i);
        end
    end
    
    % --- Nonsymmetric Hessian part
    if ~obj.ApproxAffineHessian
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
                    clear dXij
                end
                clear Bij
            end
        end
        h = (h + h')/2;   % Insure symmetric Hessian
    end
    
end
    