function [ok, model, dat] = lsSigma(ds, model, dat, opt)
% FORMAT [ok, model, dat] = lsSigma(dw, model, dat, opt)
%
% ** Required **
% dw    - Search direction
% model - Structure containing previous model data (mu, w, ww, etc.)
% dat   - Structure array containing previous individual data (f, ipsi, etc.)
% opt   - Structure containing option values (lsit, itgr, etc.)
% ** Output **
% ok    - True if a better parameter value was found
% model - Updated model
% dat   - Updated individual data
%
% Performs a line search along a direction to find a better residual 
% variance. The search direction is usually found by Gauss-Newton.

    % --- Parse inputs
    p = inputParser;
    p.FunctionName = 'lsSigma';
    p.addRequired('dw',     @checkarray);
    p.addRequired('model',  @isstruct);
    p.addRequired('dat',    @isstruct);
    p.addRequired('opt',    @isstruct);
    p.parse(ds, model, dat, opt);
    
    if opt.debug, fprintf('* lsSigma\n'); end;
    
    % --- Initial point
    sigma0 = model.sigma;
    s0     = log(sigma0);
    
    % --- Initialise line search
    armijo = 1;
    llm0   = model.llm;
    lls0   = llInvChi2((sigma0)^2, opt.sigmadf, 'fast');
    ll0    = llm0 + lls0;
    ok     = false;
    
    if opt.verbose
        printInfo('header');
        printInfo('initial', ll0, llm0, lls0);
        opt1 = opt;
        opt1.verbose = false;
    end
        
    % --- Line search
    for i=1:opt.lsit
        
        s = s0 + ds/armijo;
        model.sigma = exp(s);
%         model.sigma = sigma0 + ds/armijo;
        
        % - Update individual matching terms
        [dat, model] = pgraBatchProcess('FitResidual', dat, model, opt);
%         dat = pgraBatchProcess('Update', dat, model, opt1, ...
%             {'v', 'iphi', 'ipsi', 'pf', 'c', 'llm'}, ...
%             'clean', {'ipsi', 'iphi'});

        % - Compute log-likelihood
        llm = 0;
        for n=1:opt.N
            llm = llm + dat(n).llm;
        end
%         lls = 0;
        lls = llInvChi2((model.sigma)^2, opt.sigmadf, 'fast');
        
        ll = llm + lls;
        if opt.verbose, printInfo(armijo, ll0, llm, lls); end;
        
        % - Test improvement
        if ll <= ll0
            if opt.verbose, printInfo('failed'); end;
            armijo = armijo * 2;
        else
            ok = true;
            if opt.verbose, printInfo('success'); end;
            break
        end
    end

    % --- Clean
    if ~ok
        printInfo('end');
        model.sigma = sigma0;
        [dat, model] = pgraBatchProcess('FitResidual', dat, model, opt);
%         dat = pgraBatchProcess('Update', dat, model, opt1, ...
%             {'v', 'ipsi', 'iphi', 'pf', 'c', 'llm'}, ...
%             'clean', {'ipsi', 'iphi'});
    else
        model.llm = llm;
        model.lls = llInvChi2((model.sigma)^2, opt.sigmadf);
    end

end

function printInfo(which, oll, llm, lls)
    if ischar(which) 
        if strcmpi(which, 'header')
            fprintf([repmat('_', 1, 100) '\n']);
            fprintf('%10s | %10s | %10s = %10s + %10s   %10s  | %10s\n', 'LS Sigma', 'Armijo', 'RLL', 'LL-Match', 'LL-S', '', 'LL-Diff');
            fprintf([repmat('-', 1, 100) '\n']);
        elseif strcmpi(which, 'initial')
            fprintf('%10s | %10s | %10g = %10g + %10g   %10s\n', 'LS Sigma', 'Initial', oll, llm, lls, '');
        elseif strcmpi(which, 'failed')
            fprintf(' | Failed\n');
        elseif strcmpi(which, 'success')
            fprintf(' | Success\n');
            fprintf([repmat('_', 1, 100) '\n']);
        elseif strcmpi(which, 'end')
            fprintf('%10s | Complete failure\n', 'LS Sigma');
            fprintf([repmat('_', 1, 100) '\n']);
        end
    else
        fprintf('%10s | %10g | %10g = %10g + %10g   %10s  | %10g ', 'LS Sigma', which, llm+lls, llm, lls, '', llm+lls-oll);
    end
end