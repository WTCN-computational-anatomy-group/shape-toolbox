function ok = gaussNewtonAffine(obj)
% FORMAT ok = obj.gaussNewtonAffine()
% Performs one Gauss-Newton update of the affine transform.

    if obj.Debug, fprintf('* gaussNewtonAffine\n'); end;
    
    % --- Compute gradient and hessian
    [g, h] = obj.computeGradHessAffine();
 
    % --- Update (estimated) uncertainty on Q
    obj.hq.dim  = size(h);
    obj.hq(:)   = h(:);
    obj.statusChanged('hq');
    
    % --- Compute descent direction
    if size(obj.mu, 3) == 1
        % 2D case > Remove unused parameters
        nq = size(g);
        sub = [1 2 4 7 8 10];
        sub = sub(sub <= nq(1));
        g = g(sub);
        h = h(sub,sub);
        dq = zeros(nq);
        g = h \ g;
        dq(sub) = g;
    else
        % 3D case > easy
        dq = h\g;
    end
    clear h g
    
    % --- Line search to avoid overshoot
    ok = lineSearchQ(obj, dq);
    
%     if ok
%         obj.updateSuffStat(); % Update sufficient statistics
%     end
    
end

function ok = lineSearchQ(obj, dq, varargin)
    
    % --- Options
    opt           = struct;
    opt.MaxIt     = obj.MaxLSIt;
    opt.Verbose   = obj.Verbose;
    opt.Graphic   = obj.Graphic;
    opt           = parse_varargin(varargin, opt);
    
    if obj.Debug, fprintf('* lineSearchQ\n'); end;
    
    % --- Load data
    if ~obj.checkarray('q')
        obj.initQ();
    end
    dq = reshape(dq, [1 1 length(dq)]);
    q0 = numeric(obj.q);     % Previous reconstructed velocity
    q0 = reshape(q0, [1 1 length(q0)]);
    
    % --- Line search
    armijo = 1;
    ok = false;
    obj.logLikelihoodMatching();
    ollq = logLikelihoodPriorAffine(obj, obj.q, true);
    onll = -(obj.llm + ollq); % Only use what's improvable by Q
    if opt.Verbose
        fprintf([repmat('-', [1 94]) '\n']);
        fprintf('Q - LineSearch | Armijo  | %12s = %12s + %12s | %12s\n', 'RNLL', 'NLL-Match', 'RNLL-Prior', 'NLL-Diff');
        fprintf([repmat('-', [1 94]) '\n']);
        fprintf('Q - LineSearch | Initial | %12.6f = %12.6f + %12.6f \n', onll, -obj.llm, -ollq);
    end
    for i=1:opt.MaxIt
        
        % - Evaluate
        q = single(q0 - dq / armijo);
        A = obj.exponentiateAffine(q);
        ipsi = obj.reconstructIPsi(obj.Mmu \ A \ obj.Mf);
        [pf, pvox] = obj.pushImage(ipsi);
        llm = logLikelihoodMatching(obj, ipsi, 'warp');
        llq = logLikelihoodPriorAffine(obj, q, true);
        nll = -(llm + llq);
        
        % - Verbose/Graphic
        if opt.Verbose
            fprintf('Q - LineSearch | Try %3d | %12.6f = %12.6f + %12.6f | %12.6f ', i, nll, -llm, -llq, nll - onll);
        end
        if isfield(opt.Graphic, 'f')
            if isa(opt.Graphic.f, 'handle')
                clf(opt.Graphic.f);
            else
                opt.Graphic.f = figure;
            end
            SectionViewer(pf, 'Parent', opt.Graphic.f);
            drawnow
        end
            
        % - Check improvement
        if nll >= onll
            % Failure, try again
            if opt.Verbose
                fprintf('| Failed\n');
            end
            armijo = armijo * 2;
        else
            % Success, write data
            if opt.Verbose
                fprintf('| Success\n');
            end
            obj.q(:)        = q(:);
            obj.A(:)        = A;
            obj.ipsi.dim    = size(ipsi);
            obj.ipsi(:)     = ipsi(:);
            obj.pf.dim      = size(pf);
            obj.pf(:)       = pf(:);
            obj.pvox.dim    = size(pvox);
            obj.pvox(:)     = pvox(:);
            % Save likelihood
            obj.llm      = llm;
            obj.statusChanged('q', 'A', 'ipsi', 'pf', 'pvox', 'llm', 'llq');
            obj.logLikelihoodPriorAffine(); % < full version;
            % OK
            ok           = true;
            break;
        end
    end
    
    obj.Graphic = opt.Graphic;
    if ~ok
        if opt.Verbose
            fprintf('Q - LineSearch | Complete failure\n');
        end
    end
    
end