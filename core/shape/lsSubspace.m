function [ok, model, dat] = lsSubspace(dw, model, dat, opt, pgra)
% FORMAT [ok, model, dat] = lsSubspace(dw, model, dat, opt)
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
% Performs a line search along a direction to find a better subspace.
% The search direction is usually found by Gauss-Newton.

    % --- Parse inputs
    p = inputParser;
    p.FunctionName = 'lsSubspace';
    p.addRequired('dw',     @checkarray);
    p.addRequired('model',  @isstruct);
    p.addRequired('dat',    @isstruct);
    p.addRequired('opt',    @isstruct);
    p.parse(dw, model, dat, opt);
    
    if opt.debug, fprintf('* lsSubspace\n'); end;
    
    % --- Initial point
    w0 = model.w;
    w0.fname = fullfile(opt.directory, 'subspace_prev.nii');
    for k=1:size(model.w, 5)
        w0(:,:,:,:,k) = model.w(:,:,:,:,k);
    end
    ww0 = model.ww;
    
    % --- Initialise line search
    if isfield(model, 'armijo'),    armijo = model.armijo;
    else                            armijo = 1; end;
    llm0 = model.llm;
    llz0 = - 0.5 * trace(model.wpz(2) * model.ww * (model.Sz + model.zz));
    if opt.nz0 == 0
        llz0 = llz0 + 0.5 * opt.N * proba('LogDet', ...
            model.wpz(1) * model.Az + model.wpz(2) * model.ww);
    else
        llz0 = llz0 + model.wpz(2) * 0.5 * opt.N * proba('LogDet', model.ww);
    end
    llw0 = - 0.5 * trace(model.ww);
    ll0  = llm0 + llz0 + llw0;
    ok   = false;
    
    if opt.verbose
        printInfo('header');
        printInfo('initial', ll0, llm0, llz0, llw0);
        opt1 = opt;
        opt1.verbose = false;
    end
        
    % --- Line search
    for i=1:opt.lsit
        
        % - Explore new W
        for k=1:opt.K
            model.w(:,:,:,:,k) = w0(:,:,:,:,k) + dw(:,:,:,:,k) / armijo;
        end
        model.ww = precisionZ(model.w, opt.vs, opt.prm);
        
        % - Update individual matching terms
        dat = batchProcess('Update', dat, model, opt1, ...
            {'v', 'ipsi', 'iphi', 'pf', 'c', 'llm'}, ...
            'clean', {'ipsi', 'iphi', 'wmu'});
        
        % - Compute log-likelihood
        llm = 0;
        for n=1:opt.N
            llm = llm + dat(n).llm;
        end
        llz = - 0.5 * trace(model.wpz(2) * model.ww * (model.Sz + model.zz));
        if opt.nz0 == 0
            llz = llz + 0.5 * opt.N * proba('LogDet', ...
                model.wpz(1) * model.Az + model.wpz(2) * model.ww);
        else
            llz = llz + model.wpz(2) * 0.5 * opt.N * proba('LogDet', model.ww);
        end
        llw = - 0.5 * trace(model.ww);
        
        ll = llm + llz + llw;
        if opt.verbose, printInfo(armijo, ll0, llm, llz, llw); end;
        
        if ll <= ll0
            if opt.verbose, printInfo('failed'); end;
            armijo = armijo * 2;
        else
            if opt.verbose, printInfo('success'); end;
            if isfield(model, 'armijo')
                model.armijo = max(0.9 * armijo, 1);
            end
            ok  = true;
            break
        end
    end

    % --- Clean
    if ~ok
        printInfo('end');
        for k=1:opt.K
            model.w(:,:,:,:,k)  = w0(:,:,:,:,k);
        end
        model.ww(:) = ww0(:);
        dat = batchProcess('Update', dat, model, opt1, ...
            {'v', 'ipsi', 'iphi', 'pf', 'c', 'llm'}, ...
            'clean', {'ipsi', 'iphi'});
        model.armijo = min(1.1 * armijo, 100);
    else
        model.ll  = ll;
        model.llm = llm;
    end
    rmarray(w0);
    rmarray(ww0);
    rmarray(dw);

end

function printInfo(which, oll, llm, llz, llw)
    if ischar(which) 
        if strcmpi(which, 'header')
            fprintf([repmat('_', 1, 100) '\n']);
            fprintf('%10s | %10s | %10s = %10s + %10s + %10s  | %10s\n', 'LS PG', 'Armijo', 'RLL', 'LL-Match', 'RLL-Z', 'RLL-W', 'LL-Diff');
            fprintf([repmat('-', 1, 100) '\n']);
        elseif strcmpi(which, 'initial')
            fprintf('%10s | %10s | %10g = %10g + %10g + %10g \n', 'LS PG', 'Initial', oll, llm, llz, llw);
        elseif strcmpi(which, 'failed')
            fprintf(' | Failed\n');
        elseif strcmpi(which, 'success')
            fprintf(' | Success\n');
            fprintf([repmat('_', 1, 100) '\n']);
        elseif strcmpi(which, 'end')
            fprintf('%10s | Complete failure\n', 'LS PG');
            fprintf([repmat('_', 1, 100) '\n']);
        end
    else
        fprintf('%10s | %10g | %10g = %10g + %10g + %10g  | %10g ', 'LS PG', which, llm+llz+llw, llm, llz, llw, llm+llz+llw-oll);
    end
end