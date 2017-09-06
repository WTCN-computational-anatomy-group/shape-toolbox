function updateAffine(obj)
    
    if obj.Debug, fprintf('* updateAffine\n'); end;

    % --- Initialize Z if needed
    if ~obj.checkarray(obj.q)
        obj.initQ();
    end
    
    % --- Gauss-Newton iterations
    for i=1:obj.MaxGNIt
        ok = gaussNewtonAffine(obj);
        if ~ok,  break; end
    end
    
    if ok
        % --- Update (estimated) uncertainty on Z
        [~, h] = obj.computeGradHessAffine();
        obj.hq.dim  = size(h);
        obj.hq(:)   = h(:);
        obj.statusChanged('hq');
        clear h
    end
    
    % --- Update full log-likelihood
    obj.logLikelihood();
    if obj.Verbose, fprintf('Model log-likelihood: %d\n', obj.ll); end;
    
end