function ll = logLikelihoodLaplace(obj, hz, hr, hq)
% FORMAT ll = obj.logLikelihoodLaplace((hz), (hr), (hq))
% (hz) - Hessian of the log-likelihood w.r.t. current Z.
%        Its inverse is the covariance of Laplace approximation of 
%        p(Z | F, Mu, W, R). [default: obj.hz]
% (hr) - Hessian of the log-likelihood w.r.t. current R.
%        Its inverse is the covariance of Laplace approximation of 
%        p(R | F, Mu, W, Z). [default: obj.hr]
% (hq) - Hessian of the log-likelihood w.r.t. current Q.
%        Its inverse is the covariance of Laplace approximation of 
%        p(q | F, Mu, W, Z). [default: obj.hq]
% (ll) - Parts of the log-likelihood which encompass Laplace approximation
%        terms. [if no argout: write to obj.lll]
%
% Parts of Laplace approximation of the posterior:
% p(F|Mu) = int_{Z,R} p(F | Mu, Z, R) dZ dR
%         = (2pi)^(-K*N*Q/2) |-H(R*,Z*,q*)|^(1/2) p(F | Mu, Z*, R*, q*) p(Z*) p(R*) p(q*)
%         = (2pi)^(-K*N*Q/2) |-H(R*)|^(1/2) |-H(Z*)|^(1/2) |-H(Q*)|^(1/2) p(F | Mu, Z*, R*, q*) p(Z*) p(R*) p(q*)
% where Z*, R*, q* define a mode of p(F|Mu,Z,R,Q) and H is the hessian of
% respectively {R,Z,Q}, Z, R qnd Q around this mode.
% This function thus returns:
% ll = -0.5 * KNQ * log(2pi) + 0.5 * log(|-H(Z*)|) + 0.5 * log(|-H(R*)|) + 0.5 * log(|-H(q*)|)

    if nargout == 0 && obj.utd.lll
        ll = obj.lll;
        return
    end

    % --- Default arguments
    if nargin < 4
        if ~obj.checkarray('hq') && obj.checkarray('q')
            [~,h] = obj.computeGradHessAffine();
            obj.hq.dim = size(h);
            obj.hq(:) = h(:);
            obj.statusChanged('hq');
            clear h
        end
        hq = obj.hq;
        if nargin < 3
            if ~obj.checkarray('hr') && obj.checkarray('r')
                [~,h] = obj.computeGradHessR();
                obj.hr.dim = size(h);
                obj.hr(:) = h(:);
                obj.statusChanged('hr');
                clear h
            end
            hr = obj.hr;
            if nargin < 2
                hz = obj.hz;
                if ~obj.checkarray('hz') && obj.checkarray('z')
                    [~,h] = obj.computeGradHessZ();
                    obj.hz.dim = size(h);
                    obj.hz(:) = h(:);
                    obj.statusChanged('hz');
                    clear h
                end
            end
        end
    end
    
    ll = 0;
    
    % --- Z part
    if obj.checkarray(hz)
        K  = size(hz, 1);
        ll = ll - logDet(hz);
        useZ = 1;
    else
        K = 1;
        useZ = 0;
    end
    
    % --- R part
    if obj.checkarray(hr)
        dim = [size(hr) 1 1 1];
        N   = prod(dim(1:4));
        ll  = ll - logDetSymTensor(hr);
        useR = 1;
    else
        N = 1;
        useR = 0;
    end
    
    % --- Q part
    if obj.checkarray(hq)
        dim = [size(hq) 1 1 1];
        Q   = prod(dim(1:4));
        ll  = ll - logDet(hq);
        useQ = 1;
    else
        Q = 1;
        useQ = 0;
    end
    
    % --- Common part
    ll = ll * (-1)^(N*useR + K*useZ + Q*useQ);
    ll = ll + N^useR * K^useZ * Q^useQ * log(2 * pi);
    
    % --- Common factor
    ll = -0.5 * ll;
    
    if nargout == 0
        obj.lll     = ll;
        obj.statusChanged('lll');
    end
    
end