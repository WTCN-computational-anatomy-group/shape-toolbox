function updateR(obj)
    
    if obj.Debug, fprintf('* updateR\n'); end;

    % --- Initialize R if needed
    if ~obj.checkarray('r')
        obj.initR();
    end
    
    % --- Gauss-Newton iterations
    for i=1:obj.MaxGNIt
        ok = gaussNewtonR(obj);
        if ~ok,  break; end
    end
    
    if ok
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