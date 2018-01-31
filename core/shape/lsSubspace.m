function [ok, model, dat] = lsSubspace(dw, model, dat, opt)
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
    
    if opt.ui.debug, fprintf('* lsSubspace\n'); end
    
    % --- Initial point
    w0 = model.pg.w;
    w0.fname = fullfile(fileparts(model.pg.w.fname), 'subspace_prev.nii');
    for k=1:size(model.pg.w, 5)
        w0(:,:,:,:,k) = model.pg.w(:,:,:,:,k);
    end
    ww0 = model.pg.ww;
    
    % --- Initialise line search
    if isfield(model.pg, 'armijo'),    armijo = model.pg.armijo;
    else,                              armijo = 1; end
    llm0 = model.lb.m.val;
    llz0 = model.lb.z.val;
    llw0 = model.lb.w.val;
    ll0  = llm0 + llz0 + llw0;
    ok   = false;
    
    if opt.ui.verbose
        printInfo('header');
        printInfo('initial', ll0, llm0, llz0, llw0);
        opt1 = opt;
        opt1.ui.verbose = false;
    end
        
    % --- Line search
    for i=1:opt.iter.ls
        
        % - Explore new W
        for k=1:opt.pg.K
            model.pg.w(:,:,:,:,k) = w0(:,:,:,:,k) + dw(:,:,:,:,k) / armijo;
        end
        model.pg.ww = precisionZ(model.pg.w, opt.tpl.vs, opt.pg.prm, opt.pg.bnd);
        
        % - Update individual matching terms
        [dat, model] = pgra_batch('LB', 'Subspace', dat, model, opt1);
        
        % - Compute log-likelihood
        llm = model.lb.m.val;
        llz = model.lb.z.val;
        llw = model.lb.w.val;
        ll = llm + llz + llw;
        
        if opt.ui.verbose, printInfo(armijo, ll0, llm, llz, llw); end
        
        if ll <= ll0
            if opt.ui.verbose, printInfo('failed'); end
            armijo = armijo * 2;
        else
            if opt.ui.verbose, printInfo('success'); end
            if isfield(model.pg, 'armijo')
                model.pg.armijo = max(0.9 * armijo, 1);
            end
            ok  = true;
            break
        end
    end

    % --- Clean
    if ~ok
        printInfo('end');
        for k=1:opt.pg.K
            model.pg.w(:,:,:,:,k)  = w0(:,:,:,:,k);
        end
        model.pg.ww = ww0;
        [dat, model] = pgra_batch('LB', 'Subspace', dat, model, opt1);
        model.pg.armijo = min(1.1 * armijo, 100);
        model.pg.ok = -1;
    else
        model.pg.ok = 1;
    end
    rmarray(w0);
    rmarray(dw);

end

function printInfo(which, oll, llm, llz, llw)
    if ischar(which) 
        if strcmpi(which, 'header')
            fprintf([repmat('_', 1, 100) '\n']);
            fprintf('%10s | %10s | %10s = %10s + %10s + %10s  | %10s\n', 'LS PG', 'Armijo', 'RLL', 'LL-Match', 'RLL-Z', 'RLL-W', 'LL-Diff');
            fprintf([repmat('-', 1, 100) '\n']);
        elseif strcmpi(which, 'initial')
            fprintf('%10s | %10s | %10.4g = %10.4g + %10.4g + %10.4g \n', 'LS PG', 'Initial', oll, llm, llz, llw);
        elseif strcmpi(which, 'failed')
            fprintf(' | :(\n');
        elseif strcmpi(which, 'success')
            fprintf(' | :D\n');
            fprintf([repmat('_', 1, 100) '\n']);
        elseif strcmpi(which, 'end')
            fprintf('%10s | Complete failure\n', 'LS PG');
            fprintf([repmat('_', 1, 100) '\n']);
        end
    else
        fprintf('%10s | %10.4g | %10.4g = %10.4g + %10.4g + %10.4g  | %10.4g ', 'LS PG', which, llm+llz+llw, llm, llz, llw, llm+llz+llw-oll);
    end
end