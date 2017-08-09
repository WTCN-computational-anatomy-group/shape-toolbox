function update(obj)
    
    if obj.Debug, fprintf('* update\n'); end;

    % --- Initialize R and Z if needed
    if ~obj.checkarray('r')
        obj.initR();
    end
    if ~obj.checkarray('z')
        obj.initZ();
    end
    
    % --- Gauss-Newton iterations
    for i=1:obj.MaxGNIt
        okZ = gaussNewtonZ(obj);
        okR = gaussNewtonR(obj);
        if (~okZ) && (~okR),  break; end
    end
    
    if okZ
        % --- Update (estimated) uncertainty on Z
        [~, h] = obj.computeGradHessZ();
        obj.hz.dim  = size(h);
        obj.hz(:)   = h(:);
        obj.statusChanged('hz');
        obj.utd.hz  = true;
        clear h
    end
    if okR
        % --- Update (estimated) uncertainty on R
        [~, h] = obj.computeGradHessR();
        obj.hr.dim  = size(h);
        obj.hr(:)   = h(:);
        obj.statusChanged('hr');
        obj.utd.hr  = true;
        clear h
    end
    
    % --- Update full log-likelihood
    obj.logLikelihood();
    if obj.Verbose, fprintf('Model log-likelihood: %d\n', obj.ll); end;
end