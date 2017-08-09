function updateZ(obj)
    
    if obj.Debug, fprintf('* updateZ\n'); end;

    % --- Initialize Z if needed
    if ~obj.checkarray('z')
        obj.initZ();
    end
    
    % --- Gauss-Newton iterations
    for i=1:obj.MaxGNIt
        ok = gaussNewtonZ(obj);
        if ~ok,  break; end
    end
    
    if ok
        % --- Update (estimated) uncertainty on Z
        [~, h] = obj.computeGradHessZ();
        obj.hz.dim  = size(h);
        obj.hz(:)   = h(:);
        obj.statusChanged('hz');
        obj.utd.hz  = true;
        clear h
    end
    
    % --- Update full log-likelihood
    obj.logLikelihood();
    if obj.Verbose, fprintf('Model log-likelihood: %d\n', obj.ll); end;
    
end