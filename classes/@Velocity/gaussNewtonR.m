function ok = gaussNewtonR(obj)
% FORMAT ok = obj.gaussNewtonR()
% Performs one Gauss-Newton update of the residual field.

    if obj.Debug, fprintf('* gaussNewtonR\n'); end;
    
    % --- Compute gradient and hessian
    [g, h] = obj.computeGradHessR();
 
    % --- Update (estimated) uncertainty on R
    obj.hr.dim  = size(h);
    obj.hr(:)   = h(:);
    obj.utd.hr  = true;
    obj.statusChanged('hr');
    
    % --- Compute descent direction
    dr = spm_diffeo('fmg', single(h), single(g), [obj.VoxelSize obj.RegParam 2 2]);
    clear h g
    
    % --- Line search to avoid overshoot
    ok = lineSearchR(obj, dr);
    
%     if ok
%         obj.updateSuffStat(); % Update sufficient statistics
%     end
    
end

function ok = lineSearchR(obj, dr, varargin)
    
    % --- Options
    opt           = struct;
    opt.MaxIt     = obj.MaxLSIt;
    opt.Verbose   = obj.Verbose;
    opt.Graphic   = obj.Graphic;
    opt           = parse_varargin(varargin, opt);
    
    if obj.Debug, fprintf('* lineSearchR\n'); end;
    
    % --- Load data
    if ~obj.checkarray('r')
        obj.initR();
    end
    r0 = single(numeric(obj.r));     % Previous reconstructed velocity
    if ~obj.checkarray('v')
        obj.exponentiateVelocity();
    end
    v0 = single(numeric(obj.v));     % Previous reconstructed velocity
    
    % --- Compute small increment
    dv = obj.SigmaR * dr;
    
    % --- Line search
    armijo = 1;
    ok = false;
    obj.logLikelihoodMatching();
    ollr = logLikelihoodPriorR(obj, obj.r, true);
    onll = -(obj.llm + ollr); % Only use what's improvable by R
    if opt.Verbose
        fprintf([repmat('-', [1 94]) '\n']);
        fprintf('R - LineSearch | Armijo  | %12s = %12s + %12s | %12s\n', 'RNLL', 'NLL-Match', 'RNLL-Prior', 'NLL-Diff');
        fprintf([repmat('-', [1 94]) '\n']);
        fprintf('R - LineSearch | Initial | %12.6f = %12.6f + %12.6f \n', onll, -obj.llm, -ollr);
    end
    for i=1:opt.MaxIt
        
        % - Evaluate
        r = single(r0 - dr / armijo);
        v = single(v0 - dv / armijo);
        iphi = obj.exponentiateVelocity(v, 'iphi');
        [pf, pvox] = obj.pushImage(iphi);
        llm = logLikelihoodMatching(obj, iphi, 'warp');
        llr = logLikelihoodPriorR(obj, r, true);
        nll = -(llm + llr);
        
        % - Verbose/Graphic
        if opt.Verbose
            fprintf('R - LineSearch | Try %3d | %12.6f = %12.6f + %12.6f | %12.6f ', i, nll, -llm, -llr, nll - onll);
        end
        if isfield(opt.Graphic, 'v')
            if isa(opt.Graphic.v, 'handle')
                try
                    clf(opt.Graphic.v);
                catch
                    opt.Graphic.v = figure;
                end
            else
                opt.Graphic.v = figure;
            end
            ImageSection(v, 'Parent', opt.Graphic.v);
            drawnow
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
            obj.statusChanged('r', 'v', 'iphi', 'pf', 'pvox', 'llm', 'llr');
            obj.r.dim       = size(r);
            obj.r(:)        = r(:);
            obj.utd.r       = true;
            obj.v.dim       = size(v);
            obj.v(:)        = v(:);
            obj.utd.v       = true;
            obj.iphi.dim    = size(iphi);
            obj.iphi(:)     = iphi(:);
            obj.utd.iphi    = true;
            obj.pf.dim      = size(pf);
            obj.pf(:)       = pf(:);
            obj.utd.pf      = true;
            obj.pvox.dim    = size(pvox);
            obj.pvox(:)     = pvox(:);
            obj.utd.pvox    = true;
            % Save likelihood
            obj.llm      = llm;
            obj.utd.llm  = true;
            obj.logLikelihoodPriorR(); % < full version;
            % OK
            ok           = true;
            break;
        end
    end
    
    obj.Graphic = opt.Graphic;
    if ~ok
        if opt.Verbose
            fprintf('R - LineSearch | Complete failure\n');
        end
    end
    
end