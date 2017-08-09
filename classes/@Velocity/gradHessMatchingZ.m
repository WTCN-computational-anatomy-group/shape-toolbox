function [g, h] = gradHessMatchingZ(obj, mu, f, c, gmu, w)
% FORMAT [g, h] = obj.gradHessMatchingZ((g), (h))
% obj - A velocity object
% (mu)  - Template image ([nx ny nz nc]) [default: obj.mu]
% (f)   - Observed image pushed into template space ([nx ny nz nc])
%         [default: obj.pf]
% (c)   - Pushed voxel count ([nx ny nz nc]) [default: obj.pvox]
% (gmu) - Template spatial gradients. [default: obj.gmu]
% (w)   - Principal subspace
% (g)   - Gradient of the objective function w.r.t. (expected) latent 
%         coordinates. Incremented by this function.
% (h)   - Hessian of the objective function w.r.t. (expected) latent 
%         coordinates. Incremented by this function.
%
% Object fields [w, mu, pf, pvox, MatchingTerm] are used by this function 
% and should thus be correctly set before call.
%
% Add the matching term part to Gradient and Hessian
% This function should only be called by computeGradHess.

    % --- Default arguments
    if nargin < 6
        w = obj.w;
        if nargin < 5
            obj.computeTemplateGrad();
            gmu = obj.gmu;
            if nargin == 1
                % Lots of dependencies
                obj.reconstructTemplate();
                mu = obj.mu;
                obj.computeTemplateGrad();
                gmu = obj.gmu;
                obj.pushImage();
                f = obj.pf;
                c = obj.pvox;
            elseif nargin < 4
                obj.pushImage();
                c = obj.pvox;
                if nargin < 3
                    f = obj.pf;
                end
            end
        end
    end

    % --- Check that all arrays are ready to be used
    if ~obj.checkarray(w) || ~obj.checkarray(mu) || ...
       ~obj.checkarray(c) || ~obj.checkarray(c)  || ...
       ~obj.checkarray(gmu)
        if obj.Debug
            warning('Cannot compute gradient and hessian: missing arrays');
        end
        return
    end
    if obj.Debug, fprintf('* gradHessMatchingZ\n'); end;

    % --- Read data
    dim         = [size(w) 1 1 1];
    dim         = dim(1:5);
    dim_lattice = dim(1:3);
%     dim_vector  = dim(4); % Unused
    dim_latent  = dim(5);
    dim_classes = size(mu, 4); 
    
    % --- Allocate arrays for grad and hessian
    g = zeros([dim_latent 1]);
    if nargout > 1
        h = zeros(dim_latent);
    end

    % Loop on Z-slices to save memory
    for x3=1:dim_lattice(3)
        
        % --- Load PC in memory
        w1 = single(w(:,:,x3,:,:));
        w = zeros([dim_lattice(1:2) 1 dim_classes dim_latent], 'single');

        % --- Multiply PCs with template gradients.
        % (This way it's done once and for all)
        for k=1:dim_latent
            for l=1:dim_classes
                gmul = single(gmu(:,:,x3,l,:));
                w(:,:,:,l,k) = - w1(:,:,:,1,k) .* gmul(:,:,:,:,1) ...
                               - w1(:,:,:,2,k) .* gmul(:,:,:,:,2) ...
                               - w1(:,:,:,3,k) .* gmul(:,:,:,:,3);
                clear gmul
            end
        end
        clear w1;
        % Here, size(w) = [nx ny 1 nclasses nlatent]

        % --- Load template and pushed image
        mu = single(mu(:,:,x3,:));   % Template
        f  = single(f(:,:,x3,:));    % Pushed image
        c  = single(c(:,:,x3,:));    % Pushed voxel count

        % --- Compute grad/hess w.r.t. the complete velocity
        % (fast version that does not perform grad multiplication)
        if nargout > 1
            [g1, h1, htype] = gradHessMatchingVel(obj, mu, f, c, []);
        else
            g1 = gradHessMatchingVel(obj, mu, f, c, []);
        end
        clear mu f c

        % --- Compute grad/hess w.r.t. the latent (low-dim) coordinates
        if nargout > 1
            [g1, h1] = vel2latGradHessMatching(w, g1, h1, htype);
        else
            g1 = vel2latGradHessMatching(w, g1);
        end
        clear w
        
        % --- Increment grad/hess
        g = g + g1;
        clear g1
        if nargout > 1
            h = h + h1;
            clear h1
        end

    end % End loop on Z-slices

    % Insure symmetric hessian
    h = (h + h')/2;
    
    % Just in case
    g(~isfinite(g)) = 0;
    h(~isfinite(h)) = 0;
end

function [g, h] = vel2latGradHessMatching(w, g, h1, htype)
% FORMAT [g, h] = vel2latGradHessMatching(w, g, h, htype)
% w     - Shape subspace basis (pre-multiplied by the template gradients)
% g     - Gradient w.r.t. the full velocity
% h     - Hessian w.r.t. the full velocity
% htype - Type of the hessian approximation w.r.t. classes:
%         'diagonal' or 'symtensor'
%
% Apply the chain rule to convert Grad/Hess w.r.t. the full velocity to
% Grad/Hes w.r.t. the latent coordinates.

    % --- Dim info
    dim         = [size(w) 1 1 1];
    dim         = dim(1:5);
    dim_lattice = dim(1:3);
    dim_classes = dim(4);   % Number of classes (Pre-multiplied W(xi,k) * -Gmu(xi))
    dim_latent  = dim(5);
    
    % --- Default arguments
    do_hessian = (nargout > 1) && (nargin > 2);
    if do_hessian && nargin < 4
        % Try to guess htype
        if issame(size(h1), size(g))
            htype = 'diagonal';
        else
            htype = 'symtensor';
        end
    end
    
    % --- Gradient
    % G = w' * g1
    w = reshape(w, [], dim_latent); % Matrix form
    g = double(w' * g(:));
    
    % --- Hessian
    if do_hessian
        h = zeros(size(g, 1), 'double');
        switch htype
            case {'symtensor'}
                % H = sum_{j,k} [ w_j' * h1_{j,k} * w_k ]
                % ^ where j, k are classes/modalities
                ind = symIndices(size(h1, 4));
                w = reshape(w, [prod(dim_lattice) dim_classes dim_latent]);
                for j=1:dim_classes
                    for k=1:dim_classes
                        wk = reshape(w(:,k,:), [], dim_latent);
                        hjk = h1(:,:,:,ind(j, k));
                        % > h1_{j,k} * w_k
                        hjk = bsxfun(@times, hjk(:), wk);
                        clear wk
                        wj = reshape(w(:,j,:), [], dim_latent);
                        % > w_j' * h1_{j,k} * w_k
                        h = h + double(wj' * hjk);
                        clear hjk wj
                    end
                end
            case {'diagonal'}
                % H = sum_k [ w_k' * h1_k * w_k ]
                % ^ where k is a class/modality
                w = reshape(w, [prod(dim_lattice) dim_classes dim_latent]);
                for k=1:dim_classes
                    wk = reshape(w(:,k,:), [], dim_latent);
                    hk = h1(:,:,:,k);
                    h = h + double(wk' * bsxfun(@times, hk(:), wk));
                    clear wk hk
                end
            otherwise
                error('Unknown hessian type.');
        end
    end
    
end
