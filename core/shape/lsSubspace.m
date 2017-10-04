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
    
    if opt.debug, fprintf('* lsSubspace\n'); end;
    
    
    % --- Initial point
    w0 = model.w;
    w0.fname = fullfile(opt.directory, 'subspace_prev.nii');
    w0(:) = model.w(:);
    ww0 = model.ww;
    
    % --- Initialise line search
    if isfield(model, 'armijo'),    armijo = model.armijo;
    else                            armijo = 1; end;
    llm0  = model.llm;
    llz0 = -0.5 * trace(model.ww * (model.S + model.zz)) - 0.5 * logDet(model.ww);
    llw0 = -0.5 * trace(model.ww);
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
            {'v', 'ipsi', 'iphi', 'pf', 'c', 'llm', 'llz'}, ...
            'clean', {'ipsi', 'iphi'});
        
        % - Compute log-likelihood
        llm = 0;
        for n=1:opt.N
            llm = llm + dat(n).llm;
        end
        llz = -0.5 * model.wpz(2) * ...
                (trace(model.ww * (model.S + model.zz)) - logDet(model.ww));
        llw = -0.5 * model.wpz(1) * trace(model.ww);
        
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
        model.w(:)  = w0(:);
        model.ww(:) = ww0(:);
        dat = batchProcess('Update', dat, model, opt1, ...
            {'v', 'ipsi', 'iphi', 'pf', 'c', 'llm', 'llz'}, ...
            'clean', {'ipsi', 'iphi'});
    else
        model.ll  = ll;
        model.llm = llm;
        model.llz = llz;
        model.llw = llw;
    end
    rmarray(w0);
    rmarray(ww0);
    rmarray(dw);

end

function printInfo(which, oll, llm, llz, llw)
    if ischar(which) 
        if strcmpi(which, 'header')
            fprintf('W - LineSearch | Armijo  | %12s = %12s + %12s + %12s | %12s\n', 'RLL', 'LL-Match', 'RLL-Prior-Z', 'RLL-Prior-W', 'LL-Diff');
        elseif strcmpi(which, 'initial')
            fprintf('W - LineSearch | Initial | %12g = %12g + %12g + %12g \n', oll, llm, llz, llw);
        elseif strcmpi(which, 'failed')
            fprintf('| Failed\n');
        elseif strcmpi(which, 'success')
            fprintf('| Success\n');
        elseif strcmpi(which, 'end')
            fprintf('W - LineSearch | Complete failure\n');
        end
    else
        fprintf('W - LineSearch | %7g | %12g = %12g + %12g + %12g | %12g ', which, llm+llz+llw, llm, llz, llw, llm+llz+llw-oll);
    end
end