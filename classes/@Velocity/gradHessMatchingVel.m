function [g, h, htype] = gradHessMatchingVel(obj, mu, f, c, gmu)
% FORMAT [g, h, htype] = obj.velGradHessMatching((mu), (f), (c), (gmu))
% obj   - Velocity object
% (mu)  - Template image ([nx ny nz nc]) [default: obj.mu]
% (f)   - Observed image pushed into template space ([nx ny nz nc])
%         [default: obj.pf]
% (c)   - Pushed voxel count ([nx ny nz nc]) [default: obj.pvox]
% (gmu) - Template spatial gradients. [default: obj.gmu]
%         If provided and empty, do not use: it's deal with outside (useful
%         when computing grad/hess w.r.t. Z)
% g     - First derivatives w.r.t. full velocity ([nx ny nz nc])
% h     - Second derivatives w.r.t. full velocity
%         (diagonal approximation: [nx ny nz nc], except for the multinomial 
%         case where it is a symmetric tensor field of virtual size 
%         [nx ny nz nc nc])
% htype - Type of hessian approximation: 'diagonal' or 'symtensor'
%
% obj.MatchingTerm:
% normal      - Gaussian noise model     (F = Mu(IPhi) + N(0,s))
% laplace     - Laplace noise model      (F = Mu(IPhi) + L(0,b))
% binomial    - True/False realisation   (F = Ber(Mu(IPhi)))
% multinomial - Multiclass realisation   (F = Cat(Mu(IPhi)))
%
% Compute gradient/hessian with respect to changes in the
% complete velocity.
% Gradients actually take the form of vector fields ([nx ny nz ngrad nc]) and 
% Hessians take the form of tensor fields ([nx ny nz ngrad nc nc]), but the 
% "vector" dimension (which consists of a pointwise multiplication with 
% the template spatial gradients) is performed outside for computational 
% reasons.
    
    % --- Default arguments
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
        
    
    % --- If gradient empty, dont consider them
    check = true;
    if ~obj.checkarray(gmu)
        gmu = [];
        check = false;
    end
        
    if check
        % --- Check everything is ready to be used
        if ~obj.checkarray(mu) || ~obj.checkarray(f) || ~obj.checkarray(c)
            warning('Cannot compute gradient and hessian: missing arrays\n');
            g = [];
            h = [];
            htype = '';
            return
        end
        if obj.Debug, fprintf('* gradHessMatchingVel\n'); end;
    end
    
    % --- Compute gradient ans hessian (select case)
    switch lower(obj.MatchingTerm)
        % E = log p(f | mu, v)
        case {'normal', 'gaussian', 'l2'}
            [g, h, htype] = ghNormal(mu, f, c, obj.Normal.s, gmu);

        case {'laplace', 'l1'}
            [g, h, htype] = ghLaplace(mu, f, c, obj.Laplace.b, gmu);

        case {'bernoulli', 'binomial', 'binary'}
            [g, h, htype] = ghBernoulli(mu, f, c, gmu);

        case {'categorical', 'multinomial'}
            [g, h, htype] = ghCategorical(mu, f, c, gmu);
            
        otherwise
            error('Unknown likelihood function.');
    end
    
    % Just in case
    g(~isfinite(g)) = 0;
    h(~isfinite(h)) = 0;
end
%==========================================================================
function [g, h, htype] = ghNormal(mu, f, c, s, gmu)
% FORMAT [g, (h, htype)] = ghNormal(mu, f, c, gmu)
% mu  - Template
% f   - Observed image pushed in the template space.
% c   - Pushed voxel count.
% gmu - Template spatial gradients.
%
% Gradient & Hessian of the negative log-likelihood of the normal matching
% term w.r.t. changes in the initial velocity.
%
% If gmu is provided
% > returns g/h w.r.t. changes in v.
% > size(g) = [nx ny nz nv]
% > size(h) = [nx ny nz nv*(nv+1)/2]
%   (where nv is the number of spatial components 1/2/3)
% Else
% > Returns class-wise "residuals" that need to be multiplied with the
%   corresponding template gradient and sumed.
% > size(g) = [nx ny nz nc]
% > size(h) = [nx ny nz nc]
    
    % Gradient/Hessian of the negative log-likelihood
    % E = -log p(f | mu)
    %   = sum{k,i} 1/(2s(k)) c(xi) (mu(k)(xi) - f(k)(xi)/c(xi))^2 + Cte
    % Gd[E](xi)  = 1/s(k) (c(xi) mu(k)(xi) - f(k)(xi)) Gd[mu(k)](xi)
    % Hdl[E](xi) = 1/s(k) c(xi) Gd[mu(k)](xi) Gl[mu(k)](xi)

    nc = size(mu, 4);
    if length(s) == 1
        s = s * ones(1, nc);
    end
    c = single(numeric(c));
    
    if isempty(gmu)
    % --- Everything low-level (need to have enough memory)
        mu = single(numeric(mu));
        f  = single(numeric(f));
        s  = reshape(1./s, [1 1 1 nc]);
        g  = (bsxfun(@times, c, mu) - f);
        g  = bsxfun(@times, s, g);
        if nargout > 1
            h  = bsxfun(@times, s, c);
            htype = 'diagonal';
        end
    else
    % --- One class at a time (more i/o, less memory)
        dim  = [size(mu) 1 1];
        dlat = dim(1:3);
        nvec = size(gmu, 5);
        
        g = zeros([dlat nvec], 'single');
        if nargout > 1
            h = zeros([dlat nvec*(nvec+1)/2], 'single');
        end
        
        for k=1:nc
            mu1 = single(mu(:,:,:,k));
            f1  = single(f(:,:,:,k));
            
            % Gradient
            r1 = 1 ./ s(k) .* (c .* mu1 - f1);
            for d=1:nvec
                g(:,:,:,d) = g(:,:,:,d) + r1 .* gmu(:,:,:,k,d);
            end
            clear r1 mu1 f1
            
            % Hessian
            if nargout > 1
                ind = symIndices(size(h, 4));
                c1 = 1 / s(k) * c;
                for d1=1:nvec
                    for d2=d1:nvec
                        h(:,:,:,ind(d1,d2)) = h(:,:,:,ind(d1,d2)) ...
                            + c1 .* gmu(:,:,:,k,d1) .* gmu(:,:,:,k,d2);
                    end
                end
                clear c1
            end
        end
        htype = 'symtensor';
    end
end
%==========================================================================
function [g, h, htype] = ghLaplace(mu, f, c, b, gmu)
% FORMAT [g, (h, htype)] = ghLaplace(mu, f, c, gmu)
% mu  - Template
% f   - Observed image pushed in the template space.
% c   - Pushed voxel count.
% gmu - Template spatial gradients.
%
% Gradient & Hessian of the negative log-likelihood of the laplace 
% matching term  w.r.t. changes in the initial velocity.
%
% If gmu is provided
% > returns g/h w.r.t. changes in v.
% > size(g) = [nx ny nz nv]
% > size(h) = [nx ny nz nv*(nv+1)/2]
%   (where nv is the number of spatial components 1/2/3)
% Else
% > Returns class-wise "residuals" that need to be multiplied with the
%   corresponding template gradient and sumed.
% > size(g) = [nx ny nz nc]
% > size(h) = [nx ny nz nc]
%   (where nc is the number of classes)
    
    % Gradient/Hessian of the negative log-likelihood
    % E = -log p(f | mu)
    %   = sum{k,i} 1/b(k) c(xi) |mu(k)(xi) - f(k)(xi)/c(xi)| + Cte
    % Gd[E](xi)  = 1/b(k) sign(c(xi) mu(k)(xi) - f(k)(xi)) Gd[mu(k)](xi)
    % Hdl[E](xi) = 1/b(k) c(xi) Gd[mu(k)](xi) Gl[mu(k)](xi)
    
    nc = size(mu, 4);
    if length(b) == 1
        b = b * ones(1, nc);
    end
    c = single(numeric(c));
    
    if isempty(gmu)
    % --- Everything low-level (need to have enough memory)
        mu = single(numeric(mu));
        f = single(numeric(f));
        b = reshape(1./b, [1 1 1 nc]);
        g = bsxfun(@times, c, mu) - f;
        g = bsxfun(@times, b, sign(g));
        if nargout > 1
            h = bsxfun(@times, b, c);
            htype = 'diagonal';
        end
    else
    % --- One class at a time (more i/o, less memory)
        dim = [size(mu) 1 1];
        dlat = dim(1:3);
        nvec = size(gmu, 5);
        
        g = zeros([dlat nvec], 'single');
        if nargout > 1
            h = zeros([dlat nvec*(nvec+1)/2], 'single');
        end
        
        for k=1:nc
            mu1 = single(mu(:,:,:,k));
            f1  = single(f(:,:,:,k));
            
            % Gradient
            r1 = 1 ./ b(k) .* sign(c .* mu1 - f1);
            for d=1:nvec
                g(:,:,:,d) = g(:,:,:,d) + r1 .* gmu(:,:,:,k,d);
            end
            clear r1 mu1 f1
            
            % Hessian
            if nargout > 1
                ind = symIndices(size(h, 4));
                c1 = 1 / b(k) * c;
                for d1=1:nvec
                    for d2=d1:nvec
                        h(:,:,:,ind(d1,d2)) = h(:,:,:,ind(d1,d2)) ...
                            + c1 .* gmu(:,:,:,k,d1) .* gmu(:,:,:,k,d2);
                    end
                end
                clear c1
            end
        end
        htype = 'symtensor';
    end
end
%==========================================================================
function [g, h, htype] = ghBernoulli(mu, f, c, ga)
% FORMAT [g, (h, htype)] = ghBernoulli(mu, f, c, ga)
% mu - Reconstructed template (probability space).
% f  - Observed image pushe din the template space.
% c  - Pushed voxel count.
% ga - Tempalte spatial gradients in the log space.
%
% Gradient & Hessian of the negative log-likelihood of the bernoulli 
% matching term  w.r.t. changes in the initial velocity.
%
% If ga is provided
% > returns g/h w.r.t. changes in v.
% > size(g) = [nx ny nz nv]
% > size(h) = [nx ny nz nv*(nv+1)/2]
%   (where nv is the number of spatial components 1/2/3)
% Else
% > Returns class-wise "residuals" that need to be multiplied with the
%   corresponding template gradient.
% > size(g) = [nx ny nz]
% > size(h) = [nx ny nz]
    
    % Only one class !
    
    % Here, the input template is actually encoded in some kind of log 
    % space (that we will note a). Actual probabilities mu are 
    % mapped back with a sigmoid function:
    % mu = exp(a) / (1 + exp(a))
    
    % Gradient/Hessian of the negative log-likelihood
    % E = - log p(f | mu, v)
    %   = sum{i} c(xi) [(pf(xi)/c(xi) -1) log(1 - mu(xi)) - (pf(xi)/c(xi)) log(mu(xi))]
    %   = sum{i} [ c(xi) log(1 + exp(a(xi))) - pf(xi) a(xi) ]
    % Gd[E](xi)  = (c(xi) mu(xi) - pf(xi)) Gd[a](xi)
    % Hdl[E](xi) = c(xi) mu(xi) (1 - mu(xi)) Gd[a](xi) Gl[a](xi)
    
    % --- Load from disk (if needed)
    mu = single(numeric(mu));
    f  = single(numeric(f));
    c  = single(numeric(c));
    
    % --- Compute residuals
    g  = (c .* mu) - f;
    if nargout > 1
        h  = c .* (mu .* (1 - mu) + 1E-3);  % How to fix 1E-3 ?
        htype = 'diagonal';
    end
    
    % --- Just multiply by the spatial gradients
    if ~isempty(ga)
        % Dim info
        dim = [size(mu) 1 1];
        dlat = dim(1:3);        % Dimensions of the lattice
        nvec = size(ga, 5);     % Number of directions with a gradient
        
        % Gradient
        g1 = zeros([dlat nvec], 'single');
        [g, g1] = deal(g1, g); % swap g and g1
        for d=1:nvec
            g(:,:,:,d) = g(:,:,:,d) - g1 .* ga(:,:,:,1,d);
        end
        clear g1
        
        % Hessian
        if nargout > 1
            h1 = zeros([dlat nvec*(nvec+1)/2], 'single');
            [h, h1] = deal(h1, h); % swap h and h1
            ind = symIndices(6);
            for d1=1:nvec
                for d2=d1:nvec
                    h(:,:,:,ind(d1,d2)) = h(:,:,:,ind(d1,d2)) ...
                        + h1 .* ga(:,:,:,1,d1) .* ga(:,:,:,1,d2);
                end
            end
            clear h1
            htype = 'symtensor';
        end
    end
end
%==========================================================================
function [g, h, htype] = ghCategorical(mu, f, c, ga)
% FORMAT [g, (h, htype)] = ghCategorical(mu, f, c, ga)
% mu - Reconstructed template (probability space).
% f  - Observed image pushe din the template space.
% c  - Pushed voxel count.
% ga - Tempalte spatial gradients in the log space.
%
% Gradient & Hessian of the negative log-likelihood of the categorical 
% matching term  w.r.t. changes in the initial velocity.
%
% If ga is provided
% > returns g/h w.r.t. changes in v.
% > size(g) = [nx ny nz nv]
% > size(h) = [nx ny nz nv*(nv+1)/2]
%   (where nv is the number of spatial components 1/2/3)
% Else
% > Returns class-wise "residuals" that need to be multiplied with the
%   corresponding template gradient and sumed.
% > size(g) = [nx ny nz nc]
% > size(h) = [nx ny nz nc*(nc+1)/2]
%   (where nc is the number of classes)

    % Here, the input template is actually encoded in some kind of log 
    % space (that we will note a). Actual probabilities mu are 
    % mapped back with a softmax function:
    % mu[k] = exp(a[k]) / (sum{l} exp(a[l]))
    % Consequently, gmu actually contains Grad(a) (instead of Grad(mu))
    
    % Gradient/Hessian of the negative log-likelihood
    % E = -log p(f | mu, v)
    %   = - sum{i,k} [ c(xi) (pf(k)(xi)/c(xi)) log(mu(k)(xi)) ]
    %   = - sum{i,k} [ pf(k)(xi) log(mu(k)(xi)) ]
    % Gd[E](xi) = sum{k} (c(k)(xi) mu(k)(xi) - pf(k)(xi)) Gd[a(k)](xi)
    % Hdl[E](xi) = sum{k,j} c(xi) mu(k)(xi) (d(j,k) - mu(j)(xi)) Gd[a(k)](xi) Gl[a(j)](xi)
    
    c = single(numeric(c));
    dim = [size(mu) 1 1];
    lat = dim(1:3);
    nc = dim(4);
        
    if isempty(ga)
    % --- Everything low-level (need to have enough memory)
        % - Load from disk
        mu = single(numeric(mu));
        f  = single(numeric(f));
        
        % - Gradient
        g = bsxfun(@times, c, mu) - f;
        
        % - Hessian
        if nargout > 1
            [ind, length] = symIndices(nc, 'n');
            h = zeros([dlat length]);

            % Diagonals of the tensors
            for k=1:nc
                h(:,:,:,ind(k,k)) = c .* (  mu(:,:,:,k) - mu(:,:,:,k).^2 );
            end

            % Upper/Lower part
            for l=1:nc
                for m=(l+1):nc
                    h(:,:,:,ind(l,m)) = - c .* mu(:,:,:,l) .* mu(:,:,:,m);
                end
            end
            htype = 'symtensor';
        end
        
    else
    % --- One class at a time (more i/o, less memory)
        g = zeros([lat 3], 'single');
        if nargout > 1
            h = zeros([lat 6], 'single');
            ind = symIndices(6);
        end
        
        for k=1:nc
            % - Load data from disk
            mu1 = single(mu(:,:,:,k));
            f1 = single(f(:,:,:,k));
            
            % - Gradient
            r1 = c.* mu1 - f1;
            clear f1
            for d=1:3
                ga1 = single(ga(:,:,:,k,d));
                g(:,:,:,d) = g(:,:,:,d) + r1 .* ga1;
                clear gmu1
            end
            
            % - Hessian
            if nargout > 1
                for kk=1:nc
                    if k == kk
                        % Diagonal
                        c1 = c .* (mu1 - mu1.^2);
                    else
                        % Off-diagonal
                        mu2 = single(mu(:,:,:,kk));
                        c1 = - c .* mu1 .* mu2;
                        clear mu2
                    end
                    % Multiply with spatial gradients
                    for d1=1:3
                        ga1 = single(ga(:,:,:,k,d1));
                        for d2=d1:3
                            ga2 = single(ga(:,:,:,k,d2));
                            h(:,:,:,ind(d1,d2)) = h(:,:,:,ind(d1,d2)) + c1 .* ga1 .* ga2;
                            clear gmu2
                        end
                        clear gmu1
                    end
                    clear c1
                end
                clear mu1 
            end
        end
        htype = 'symtensor';
    end
        
end