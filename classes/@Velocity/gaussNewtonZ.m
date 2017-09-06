function ok = gaussNewtonZ(obj)
% FORMAT ok = objgaussNewtonZ()
% Performs one Gauss-Newton update of the latent coordinates.

    if obj.Debug, fprintf('* gaussNewtonZ\n'); end;
    
    % --- Compute gradient and hessian
    [g, h] = obj.computeGradHessZ();
    
    % --- Update (estimated) uncertainty on Z
    obj.hz.dim  = size(h);
    obj.hz(:)   = h(:);
    obj.utd.hz  = true;
    obj.statusChanged('hz');
    
    % --- Compute descent direction
    dz = h\g;
    clear h g
 
    % --- Line search to avoid overshoot
    ok = lineSearchZ(obj, dz);
    
%     if ok
%         obj.updateSuffStat(); % Update sufficient statistics
%     end
    
end

function ok = lineSearchZ(obj, dz, varargin)
    
    % --- Options
    opt           = struct;
    opt.MaxIt     = obj.MaxLSIt;
    opt.Verbose   = obj.Verbose;
    opt.Graphic   = obj.Graphic;
    opt           = parse_varargin(varargin, opt);
    
    if obj.Debug, fprintf('* lineSearchZ\n'); end;
    
    % --- Load data
    if ~obj.checkarray('z')
        obj.initZ();
    end
    z0 = single(numeric(obj.z));     % Previous latent coordinates
    if ~obj.checkarray('v')
        obj.exponentiateVelocity();
    end
    v0 = single(numeric(obj.v));     % Previous reconstructed velocity
    
    % --- Compute small increment
    dv = lat2vel(dz, obj.w);
    
    % --- Line search
    armijo = 1;
    ok = false;
    obj.logLikelihoodMatching();
    ollz = logLikelihoodPriorZ(obj, obj.z, obj.regz, true); % < fast version
    onll = -(obj.llm + ollz); % Only use what's improvable by Z
    if opt.Verbose
        fprintf('-------------------------------------------------------\n');
        fprintf('Z - LineSearch | Armijo  | %12s = %12s + %12s | %12s\n', 'RNLL', 'NLL-Match', 'RNLL-Prior', 'NLL-Diff');
        fprintf('-------------------------------------------------------\n');
        fprintf('Z - LineSearch | Initial | %12.8f = %12.8f + %12.8f \n', onll, -obj.llm, -ollz);
    end
    for i=1:opt.MaxIt
        
        % - Evaluate
        z = single(z0 - dz / armijo);
        v = single(v0 - dv / armijo);
        iphi = obj.exponentiateVelocity(v, 'iphi');
        [pf, pvox] = obj.pushImage(iphi);
        llm = logLikelihoodMatching(obj, iphi, 'warp');
        llz = logLikelihoodPriorZ(obj, z, obj.regz, true); % < fast version
        nll = -(llm + llz);
        
        % - Verbose/Graphic
        if opt.Verbose
            fprintf('Z - LineSearch | Try %3d | %12.8f = %12.8f + %12.8f | %12.8f ', armijo, nll, -llm, -llz, nll - onll);
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
            obj.z.dim       = size(z);
            obj.z(:)        = z(:);
            obj.v.dim       = size(v);
            obj.v(:)        = v(:);
            obj.iphi.dim    = size(iphi);
            obj.iphi(:)     = iphi(:);
            obj.pf.dim      = size(pf);
            obj.pf(:)       = pf(:);
            obj.pvox.dim    = size(pvox);
            obj.pvox(:)     = pvox(:);
            % Save likelihood
            obj.llm      = llm;
            obj.utd.llm  = true;
            obj.statusChanged('z', 'v', 'iphi', 'pf', 'pvox', 'llm');
            obj.logLikelihoodPriorZ(); % < full version;
            % OK
            ok           = true;
            break;
        end
    end
    
    obj.Graphic = opt.Graphic;
    if ~ok
        if opt.Verbose
            fprintf('Z - LineSearch | Complete failure\n');
        end
    end
    
end