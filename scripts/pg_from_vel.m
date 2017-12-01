function [model, dat] = pg_from_vel(opt, dat, model, cont)
% _________________________________________________________________________
%
% FORMAT [model, dat] = pg_from_vel(opt, dat, model)
%
% Learn a principal subspace from velocity fields.
%
% The input option structure (opt) should at least contain the field:
% fnames.dat.v - observed velocity fields (as a list of filenames)
%
% Alternatively, an input data structure array (dat), of size the number
% of images, can be provided with at least the field:
% v        - observed velocity fields (as a list of array or file_array)
%
% All velocity fields should have identical dimensions and voxel sizes.
%
% The following parameters can be overriden by specifying them in the
% input option (opt) structure:
% ** MODEL **
% K            - Number of principal geodesics [32]
% vs           - Lattice voxel size [auto]
% prm          - Parameters of the geodesic differential operator
%                [1e-4 1e-3 0.2 0.05 0.2]
% Az0          - Prior expected value of the latent precision matrix [eye(K)] 
% nz0          - Prior degrees of freedom of the latent precision matrix [K]
% lam0         - Prior expected value of the anatomical noise precision [17]
% nl0          - Prior degrees of freedom of the anatomical noise precision [10]
% ** PROCESSING **
% emit         - Maximum number of EM iterations [100]
% itgr         - Number of integration steps for geodesic shooting [auto]
% verbose      - Talk during line search [true]
% debug        - Further debuging talk [false]
% par          - Parallelise processing 0/n/inf [inf]
% batch        - Batch size for parallelisation [auto].
% ** I/O **
% directory    - Directory where to store result arrays ['.']
% fnames.result- Filename for the result environment saved after each EM
%                iteration ['pg_result.mat']
% fnames.model - Structure of filenames for all temporary arrays
%                (mu, gmu, (a), w, dw, g, h)
% fnames.dat   - Structure of filenames for all temporary arrays
%                (f, v, ipsi, iphi, pf, c, wmu, r)
% ondisk.model - Structure of logical for temporary array [default_ondisk]
% ondisk.dat   - "      "       "       "       "       "
% _________________________________________________________________________
%
% FORMAT [model, dat, opt] = pg_from_vel(opt, dat, model, 'continue')
%
% The returned structures (or the saved environment which also contains
% them) can be used as input to start optimising from a previous state.
% _________________________________________________________________________

% _________________________________________________________________________
%
% Graphical model
% ---------------
%                          [L] --> (w)
%                                   |
%                                   v
%  [Az0, nz0] --> (Az) --> (z) --> {v} <-- (lam) <-- [lam0, nl0]
%                                   ^
%                                   |
%                                  [L]
% _________________________________________________________________________

    % ---------------------------------------------------------------------
    %    Default parameters
    % ---------------------------------------------------------------------
    if nargin < 4
        cont = '';
        if nargin < 3
            model = struct;
            if nargin < 2
                dat = struct;
                if nargin < 1
                    opt = struct;
                end
            end
        end
    end
    cont = strcmpi(cont, 'continue');
    
    [opt, dat]        = pg_from_vel_input(opt, dat);
    opt               = pg_from_vel_default(opt);
    [opt, dat, model] = pg_from_vel_data(opt, dat, model);
    
    opt.lb.track = true;
    
    
    if opt.N < opt.K
        warning(['It makes no sense to compute more principal components ' ...
                 'than subjects. Setting K = N-1 = %d.'], opt.N - 1);
        opt.K   = opt.N - 1;
        opt.Az0 = opt.Az0(1:opt.K,1:opt.K);
    end
    
    % ---------------------------------------------------------------------
    %    Initialise all variables
    % ---------------------------------------------------------------------
    if ~cont
        if opt.verbose
            fprintf(['%10s | %10s | ' repmat('=',1,50) ' |\n'], 'EM', 'Init');
        end
        [dat, model] = initAll(dat, model, opt);
    end
    
    % ---------------------------------------------------------------------
    %    EM iterations
    % ---------------------------------------------------------------------
    for emit = 1:opt.emit
    
        if opt.verbose
            fprintf(['%10s | %10d | ' repmat('=',1,50) ' |\n'], 'EM', emit);
        end
        
        % -----------------------------------------------------------------
        %    Principal subspace
        % -----------------------------------------------------------------
        
        if opt.verbose, fprintf('%10s | %10s ', 'PG', ''); tic; end;
        M = model.Sz + model.zz + eye(opt.K) / model.lam;
        M = loadDiag(M);
        P = model.Z' / M;
        % Lots of i/o but no explosion of memory
        for k=1:opt.K
            model.w(:,:,:,:,k) = 0;
        end
        for n=1:opt.N
            v1 = numeric(dat(n).v);
            for k=1:opt.K
                model.w(:,:,:,:,k) = model.w(:,:,:,:,k) + v1 * P(n,k);
            end
        end
        clear M P
        model.ww = precisionZ(model.w, opt.vs, opt.prm);
        if opt.verbose, fprintf('| %6gs\n', toc); end;
        
        % -----------
        % Lower bound
        if opt.lb.track
            [dat, model]    = batchLLV(dat, model, opt);
            model.lb.v.list = [model.lb.v.list model.lb.v.val];
            model.lb.v.it   = [model.lb.v.it   emit];
            model.lb.w.val  = llPriorSubspace(model.w, model.ww, opt.vs, opt.prm);
            model.lb.w.list = [model.lb.w.list model.lb.w.val];
            model.lb.w.it   = [model.lb.w.it   emit];
            model           = updateLowerBound(model);
        end
        plotAll(model, opt);
        % -----------
            
        % -----------------------------------------------------------------
        %    Latent coordinates
        % -----------------------------------------------------------------
        
        % Coordinates (variational update)
        % --------------------------------
        [dat, model] = batchLatent(dat, model, opt);

        % -----------
        % Lower bound
        if opt.lb.track
            [dat, model]    = batchLLV(dat, model, opt);
            model.lb.v.list = [model.lb.v.list model.lb.v.val];
            model.lb.v.it   = [model.lb.v.it   emit + 1/3];
            model.lb.z.list = [model.lb.z.list model.lb.z.val];
            model.lb.z.it   = [model.lb.z.it   emit];
            model           = updateLowerBound(model);
        end
        plotAll(model, opt);
        % -----------
        
%         % Orthogonalisation matrix
%         % ------------------------
%         if opt.verbose, fprintf('%10s | %10s ', 'Ortho', ''); tic; end;
%         [U, iU] = orthogonalisationMatrix(model.zz + model.Sz, model.ww);
%         if opt.verbose, fprintf('| %6gs\n', toc); end;
%         
%         % Rescaling
%         % ---------
%         if opt.verbose, fprintf('%10s | %10s ', 'Rescale', ''); tic; end;
%         zz = U*model.zz*U';
%         Sz = U*model.Sz*U';
%         ww = iU'*model.ww*iU;
%         [Q, iQ] = gnScalePG_vel(Sz, zz, ww, opt.nz0, opt.N);
% %         Q  = eye(opt.K);
% %         iQ = Q;
%         if opt.verbose, fprintf('| %6gs\n', toc); end;
%         
%         % Apply full transform
%         % --------------------
%         Q  = Q*U;
%         iQ = iU*iQ;
%         [dat, model] = batchLatentRotate(dat, model, opt, Q);
%         for z=1:size(model.w, 3)
%             w1 = reshape(model.w(:,:,z,:,:), [], opt.K);
%             w1 = w1 * iQ;
%             w1 = reshape(w1, [opt.lat(1:2) 1 3 opt.K]);
%             model.w(:,:,z,:,:) = w1;
%         end
%         model.ww = iQ' * model.ww * iQ;
        
        % Latent precision
        % ----------------
        if opt.verbose, fprintf('%10s | %10s ', 'Lat Prec', ''); tic; end;
        model.Az = (opt.nz0 + opt.N) * inv(opt.nz0 * inv(opt.Az0) + model.Sz + model.zz);
        if opt.verbose, fprintf('| %6gs\n', toc); end;
        
        % -----------
        % Lower bound
        if opt.lb.track
            model.lb.v.list = [model.lb.v.list model.lb.v.val];
            model.lb.v.it   = [model.lb.v.it   emit + 0.5];
            model.lb.z.list = [model.lb.z.list model.lb.z.val];
            model.lb.z.it   = [model.lb.z.it   emit];
            model.lb.a.val  = -klWishart(opt.nz0 + opt.N, model.Az, ...
                                         opt.nz0, opt.Az0, ...
                                         'expectation');
            model.lb.a.list = [model.lb.a.list model.lb.a.val];
            model.lb.a.it   = [model.lb.a.it   emit];
            model           = updateLowerBound(model);
        end
        plotAll(model, opt);
        % -----------
        
        if true
        % -----------------------------------------------------------------
        %    Residual precision
        % -----------------------------------------------------------------
        
        if opt.verbose, fprintf('%10s ', 'Res Prec'); tic; end;
        model.lam = opt.nl0/opt.lam0 + (trace(model.Sz*model.ww) + model.cumres)/(prod(opt.lat)*3);
        model.lam = (opt.nl0 + opt.N)/model.lam;
        model.lami = [model.lami emit];
        model.lams = [model.lams model.lam];
        if opt.verbose, fprintf('| %10.3f | %6gs\n', model.lam, toc); end;
        
        % -----------
        % Lower bound
        if opt.lb.track
            [dat, model] = batchLLV(dat, model, opt);
            model.lb.v.list = [model.lb.v.list model.lb.v.val];
            model.lb.v.it   = [model.lb.v.it   emit + 2/3];
            model.lb.l.val  = -klGamma((opt.nl0 + opt.N), model.lam, ...
                                       opt.nl0, opt.lam0, ...
                                       'precision', prod(opt.lat)*3);
            model.lb.l.list = [model.lb.l.list model.lb.l.val];
            model.lb.l.it   = [model.lb.l.it   emit];
            model           = updateLowerBound(model, true); % COMPUTE GAIN
        end
        plotAll(model, opt);
        % -----------
        end
        
        % -----------------------------------------------------------------
        %    Save current state
        % -----------------------------------------------------------------
        if ~isempty(opt.fnames.result)
            save(fullfile(opt.directory, opt.fnames.result), ...
                 'model', 'dat', 'opt');
        end
        
        if opt.lb.track && model.lb.lb.gain < 1e-7
            fprintf('Converged :D\n');
            break;
        end
    end
end

% =========================================================================

function model = updateLowerBound(model, last)
% FORMAT model = updateLowerBound(model, (last))
% model - Model structure
% last  - Last call before new EM loop ? [false]
%
% Compute the lower bound from its subparts
% + Compute gain against the previous VEM iteration if [last == true]
% + Some stuff to have nice plots where x values are actual EM iterations.

    if nargin < 2
        last = false;
    end

    % Initialise
    % ----------
    if ~isfield(model.lb, 'lb')
        model.lb.lb.val     = -inf;
        model.lb.lb.list    = []; % LB values from previous iterations
        model.lb.lb.it      = []; % x-axis values from previous iterations
        model.lb.lb.curlist = []; % LB values from current iteration
        model.lb.lb.curit   = []; % x-axis values from current iteration
    end
    
    % Compute current value
    % ---------------------
    preval = model.lb.lb.val;
    vars = fieldnames(model.lb);
    model.lb.lb.val = 0;
    for i=1:numel(vars)
        var = vars{i};
        if ~any(strcmpi(var, {'track', 'lb'}))
            model.lb.lb.val = model.lb.lb.val + model.lb.(var).val;
        end
    end
    model.lb.lb.curlist = [model.lb.lb.curlist model.lb.lb.val];
    diff = model.lb.lb.val - preval;

    % Update iterations (for x axis)
    % ------------------------------
    if isempty(model.lb.lb.it)
        it = 0;
    else
        it = floor(model.lb.lb.it(end)) + 1;
    end
    N  = length(model.lb.lb.curlist);
    model.lb.lb.curit = it + (0:(N-1))/N;
    
    % Compute gain
    % ------------
    if last
        if isempty(model.lb.lb.list)
            model.lb.lb.gain = inf;
        else
            model.lb.lb.gain = (model.lb.lb.val - model.lb.lb.list(end))/abs(model.lb.lb.list(end));
        end
        model.lb.lb.list    = [model.lb.lb.list model.lb.lb.curlist];
        model.lb.lb.it      = [model.lb.lb.it   model.lb.lb.curit];
        model.lb.lb.curlist = [];
        model.lb.lb.curit   = [];
    end
    
    % Print stuff
    % -----------
    fprintf('%10s | ', 'LB');
    if diff > 0
        fprintf('%10s | ', '(+)');
    elseif diff < 0
        fprintf('%10s | ', '(-)');
    else
        fprintf('%10s | ', '(=)');
    end
    fprintf(' %6g', model.lb.lb.val);
    if last
        fprintf([repmat(' ', 1, 37) ' | %6e'], model.lb.lb.gain);
    end
    fprintf('\n')
end

% =========================================================================

function plotAll(model, opt)
% Plot PG + lower bound stuff
% This function is highly specific to this particular model. Not sure I can
% come up with a generic "plotModel" function, even though it could be
% nice.
    
    if opt.verbose
        nw = 3;
        nh = 4;
        i  = 1;
        colors = ['b', 'g', 'r', 'c', 'm', 'k'];
        
        % PG 1/2/3
        [X,Y] = ndgrid(1:opt.lat(1), 1:opt.lat(2));
        Z = ceil(opt.lat(3)/2);
        for k=1:3
            U = reshape(model.w(:,:,Z,1,k), [opt.lat(1) opt.lat(2)]);
            V = reshape(model.w(:,:,Z,1,k), [opt.lat(1) opt.lat(2)]);
            subplot(nh,nw,i)
            quiver(X,Y,U,V)
            axis off
            title(sprintf('PG %d', k))
            i = i + 1;
        end
        
        % Lower bound
        subplot(nh,nw,i)
        plot([model.lb.lb.it model.lb.lb.curit], ...
             [model.lb.lb.list model.lb.lb.curlist], ...
             colors(mod(i, length(colors))+1))
        title('Lower bound')
        i = i + 1;
        % LL velocity
        subplot(nh,nw,i)
        plot(model.lb.v.it, model.lb.v.list, ...
             colors(mod(i, length(colors))+1))
        title('LL velocity')
        i = i + 1;
        % LL PG
        subplot(nh,nw,i)
        plot(model.lb.w.it, model.lb.w.list, ...
             colors(mod(i, length(colors))+1))
        title('LL subspace')
        i = i + 1;
        % KL latent
        subplot(nh,nw,i)
        plot(model.lb.z.it, model.lb.z.list, ...
             colors(mod(i, length(colors))+1))
        title('-KL latent')
        i = i + 1;
        % KL latent precision
        subplot(nh,nw,i)
        plot(model.lb.a.it, model.lb.a.list, ...
             colors(mod(i, length(colors))+1))
        title('-KL latent precision')
        i = i + 1;
        % KL residual precision
        subplot(nh,nw,i)
        plot(model.lb.l.it, model.lb.l.list, ...
             colors(mod(i, length(colors))+1))
        title('-KL residual precision')
        i = i + 1;
        
        % Precision matrix
        subplot(nh,nw,i)
        imagesc(model.Az), colorbar;
        title('Precision matrix')
        i = i + 1;
        % Precision WW
        subplot(nh,nw,i)
        imagesc(model.ww), colorbar;
        title('W''LW')
        i = i + 1;
        
        % Lambda
        subplot(nh,nw,i)
        plot(model.lami, model.lams, ...
             colors(mod(i, length(colors))+1));
        title('Lambda')
        
        drawnow
    end

end

% =========================================================================

function [dat, model] = initAll(dat, model, opt)
% Initialise all variables (that need it) 
% + lower bound stuff

    % Momentum
    % --------
    dat = batchMomentum(dat, opt);

    % Principal subspace
    % ------------------
    if ~isfield(opt, 'init') || ~isfield(opt.init, 'w') || opt.init.w
        model.w.dim = [opt.lat 3 opt.K];
        model.w(:)  = 0;
        model.ww    = zeros(opt.K);
    end
    if opt.lb.track
        model.lb.w.val  = llPriorSubspace(model.w, model.ww, opt.vs, opt.prm);
        model.lb.w.list = model.lb.w.val;
        model.lb.w.it   = 0;
        model.lb.w.type = 'll';
        model.lb.v.name = 'Subspace';
    end
    
    % Latent precision
    % ----------------
    model.Az = opt.Az0;
    if opt.lb.track
        model.lb.a.val  = 0;
        model.lb.a.list = 0;
        model.lb.a.it   = 0;
        model.lb.a.type = 'kl';
        model.lb.v.name = 'Latent precision';
    end
    
    % Residual precision
    % ------------------
    model.lam  = opt.lam0;
    model.lams = opt.lam0;
    model.lami = 0;
    if opt.lb.track
        model.lb.l.val  = 0;
        model.lb.l.list = 0;
        model.lb.l.it   = 0;
        model.lb.l.type = 'kl';
        model.lb.v.name = 'Residual precision';
    end
    
    % Latent coordinates
    % ------------------
    if ~isfield(opt, 'init') || ~isfield(opt.init, 'z') || opt.init.z
        [dat, model] = batchLatentInit(dat, model, opt);
    end
    if opt.lb.track
        model.lb.z.list = model.lb.z.val;
        model.lb.z.it   = 0;
        model.lb.z.type = 'kl';
        model.lb.v.name = 'Latent';
    end

    % Velocity
    % --------
    if opt.lb.track
        [dat, model]    = batchLLV(dat, model, opt);
        model.lb.v.list = model.lb.v.val;
        model.lb.v.it   = 0;
        model.lb.v.type = 'll';
        model.lb.v.name = 'Velocity';
    end
    
    % Lower Bound
    % -----------
    model = updateLowerBound(model, true);
    plotAll(model, opt);
end

% =========================================================================

function dat = oneStepMomentum(dat, opt)
% Compute momentum from velocity once and for al

    mom = spm_diffeo('vel2mom', single(numeric(dat.v)), double([opt.vs, opt.prm]));
    dat.m.dim = size(dat.v);
    dat.m(:,:,:,:) = mom(:,:,:,:);
    clear mom
end

function dat = batchMomentum(dat, opt)

    % --- Detect parallelisation scheme
    if opt.par > 0
        batch = opt.batch;
    else
        batch = 1;
    end
    
    % --- Batch processing
    if opt.verbose, before = plotBatchBegin('Momentum'); end;
    for i=1:ceil(opt.N/batch)
        n1 = (i-1)*batch + 1;
        ne = min(opt.N, i*batch);

        if opt.verbose, before = plotBatch(i, batch, opt.N, 50, before); end;
    
        % Compute subjects log-likelihood
        % -------------------------------
        parfor (n=n1:ne, double(opt.par))
            dat(n) = oneStepMomentum(dat(n), opt);
        end
        
    end
    if opt.verbose, plotBatchEnd; end;

end

% =========================================================================

function dat = oneStepLLV(dat, model, opt)
% Compute Expected log-likelihood of the data
% + store (v-Wz)'L(v-Wz) (useful for later)
    res = numeric(dat.v) - reconstructVelocity('latent', dat.z, 'subspace', model.w);
    mom = spm_diffeo('vel2mom', single(res), double([opt.vs, opt.prm]));
    dat.res = res(:)' * mom(:);
    clear res mom
    
    M = prod(opt.lat) * 3;
    dat.ll = - M * log(2*pi) ...
             + M * (log(model.lam) + ...
                    psi(0.5*(opt.nl0 + opt.N)*M) - ...
                    log(0.5*(opt.nl0 + opt.N)*M)) ...
             + proba('LogDetDiffeo', opt.lat, opt.vs, opt.prm) ...
             - model.lam * dat.res ...
             - model.lam * trace(dat.Sz * model.ww);
    dat.ll = dat.ll * 0.5;
end

function [dat, model] = batchLLV(dat, model, opt)

    % --- Detect parallelisation scheme
    if opt.par > 0
        batch = opt.batch;
    else
        batch = 1;
    end
    
    % Init accu
    % ---------
    model.cumres   = 0;
    model.lb.v.val = 0;
    
    % --- Batch processing
    if opt.verbose, before = plotBatchBegin('LL Vel'); end;
    for i=1:ceil(opt.N/batch)
        n1 = (i-1)*batch + 1;
        ne = min(opt.N, i*batch);

        if opt.verbose, before = plotBatch(i, batch, opt.N, 50, before); end;
    
        % Compute subjects log-likelihood
        % -------------------------------
        parfor (n=n1:ne, double(opt.par))
            dat(n) = oneStepLLV(dat(n), model, opt);
        end
        
        for n=n1:ne
            
            % Add individual contributions
            % ----------------------------
            model.lb.v.val = model.lb.v.val + dat(n).ll;
            model.cumres   = model.cumres   + dat(n).res;
        end
        
    end
    if opt.verbose, plotBatchEnd; end;

end

% =========================================================================

function dat = oneStepLatent(dat, model, opt)
% Update mean and covariance of the latent coordinates
% + KL divergence

    prec = model.lam * model.ww + model.Az;
    
    mom = numeric(dat.m);
    for k=1:opt.K
        w1 = single(model.w(:,:,:,:,k));
        dat.z(k) = w1(:)' * mom(:);
    end
    clear w1 mom
    dat.z = model.lam * (prec \ dat.z);
    dat.Sz = inv(prec);
    
    dat.klz = trace(prec\model.Az) ...
              - proba('LogDet', model.Az) ...
              + proba('LogDet', prec) ...
              + dat.z' * model.Az * dat.z ...
              - opt.K;
    dat.klz = 0.5 * dat.klz;
    
end

function [dat, model] = batchLatent(dat, model, opt)

    % --- Detect parallelisation scheme
    if opt.par > 0
        batch = opt.batch;
    else
        batch = 1;
    end
    
    % Init accu
    % ---------
    model.lb.z.val = 0;
    model.Sz = zeros(opt.K);
    model.zz = zeros(opt.K);
    model.z  = zeros(opt.K, 1);
    model.Z  = zeros(opt.K, opt.N);
    
    % --- Batch processing
    if opt.verbose, before = plotBatchBegin('Latent'); end;
    for i=1:ceil(opt.N/batch)
        n1 = (i-1)*batch + 1;
        ne = min(opt.N, i*batch);

        if opt.verbose, before = plotBatch(i, batch, opt.N, 50, before); end;
    
        % Compute subjects log-likelihood
        % -------------------------------
        parfor (n=n1:ne, double(opt.par))
            dat(n) = oneStepLatent(dat(n), model, opt);
        end
        
        for n=n1:ne
            
            % Add individual contributions
            % ----------------------------
            model.lb.z.val = model.lb.z.val - dat(n).klz;
            model.Sz       = model.Sz       + dat(n).Sz;
            model.zz       = model.zz       + dat(n).z * dat(n).z';
            model.z        = model.z        + dat(n).z;
            model.Z(:,n)   = dat(n).z;
        end
        
    end
    if opt.verbose, plotBatchEnd; end;

end


% =========================================================================

function dat = oneStepLatentRotate(dat, model, opt, T, iT)
% Update mean and covariance of the latent coordinates
% + KL divergence

    dat.z = T * dat.z;
    dat.Sz = T * dat.Sz * T';

    prec = model.lam * model.ww + model.Az;
    dat.klz = trace(prec\model.Az) ...
              - proba('LogDet', model.Az) ...
              + proba('LogDet', prec) ...
              + dat.z' * model.Az * dat.z ...
              - opt.K;
    dat.klz = 0.5 * dat.klz;
    
end

function [dat, model] = batchLatentRotate(dat, model, opt, T, iT)

    if nargin < 5
        iT = inv(T);
    end

    % --- Detect parallelisation scheme
    if opt.par > 0
        batch = opt.batch;
    else
        batch = 1;
    end
    
    % Init accu
    % ---------
    model.lb.z.val = 0;
    model.Sz = zeros(opt.K);
    model.zz = zeros(opt.K);
    model.z  = zeros(opt.K, 1);
    model.Z  = zeros(opt.K, opt.N);
    
    % --- Batch processing
    if opt.verbose, before = plotBatchBegin('Rotate'); end;
    for i=1:ceil(opt.N/batch)
        n1 = (i-1)*batch + 1;
        ne = min(opt.N, i*batch);

        if opt.verbose, before = plotBatch(i, batch, opt.N, 50, before); end;
    
        % Compute subjects log-likelihood
        % -------------------------------
        parfor (n=n1:ne, double(opt.par))
            dat(n) = oneStepLatentRotate(dat(n), model, opt, T, iT);
        end
        
        for n=n1:ne
            
            % Add individual contributions
            % ----------------------------
            model.lb.z.val = model.lb.z.val - dat(n).klz;
            model.Sz       = model.Sz       + dat(n).Sz;
            model.zz       = model.zz       + dat(n).z * dat(n).z';
            model.z        = model.z        + dat(n).z;
            model.Z(:,n)   = dat(n).z;
        end
        
    end
    if opt.verbose, plotBatchEnd; end;

end

% =========================================================================

function dat = oneStepLatentKL(dat, model, opt)
% Update KL divergence of the latent coordinates
    
    dat.klz = trace(dat.Sz*model.Az) ...
              - proba('LogDet', model.Az) ...
              - proba('LogDet', dat.Sz) ...
              + dat.z' * model.Az * dat.z ...
              - opt.K;
    dat.klz = 0.5 * dat.klz;
    % (We actually store -KL)
    
end

function [dat, model] = batchLatentKL(dat, model, opt)

    % --- Detect parallelisation scheme
    if opt.par > 0
        batch = opt.batch;
    else
        batch = 1;
    end
    
    % Init accu
    % ---------
    model.lb.z.val = 0;
    
    % --- Batch processing
    if opt.verbose, before = plotBatchBegin('KL Latent'); end;
    for i=1:ceil(opt.N/batch)
        n1 = (i-1)*batch + 1;
        ne = min(opt.N, i*batch);

        if opt.verbose, before = plotBatch(i, batch, opt.N, 50, before); end;
    
        % Compute subjects log-likelihood
        % -------------------------------
        parfor (n=n1:ne, double(opt.par))
            dat(n) = oneStepLatentKL(dat(n), model, opt);
        end
        
        for n=n1:ne
            
            % Add individual contributions
            % ----------------------------
            model.lb.z.val = model.lb.z.val - dat(n).klz;
        end
        
    end
    if opt.verbose, plotBatchEnd; end;

end


% =========================================================================

function dat = oneStepLatentInit(dat, model, opt)
% Update mean and covariance of the latent coordinates
% + KL divergence

    prec   = model.lam * model.ww + model.Az;
    dat.z  = mvnrnd(zeros(opt.K,1), model.Az);
    dat.z  = dat.z(:);
    dat.Sz = inv(prec);
    
    dat.klz = trace(prec\model.Az) ...
              - proba('LogDet', model.Az) ...
              + proba('LogDet', prec) ...
              + dat.z' * model.Az * dat.z ...
              - opt.K;
    dat.klz = -0.5 * dat.klz;
    % (We actually store -KL)
    
end

function [dat, model] = batchLatentInit(dat, model, opt)

    % --- Detect parallelisation scheme
    if opt.par > 0
        batch = opt.batch;
    else
        batch = 1;
    end
    
    % Init accu
    % ---------
    model.lb.z.val = 0;
    model.Sz = zeros(opt.K);
    model.zz = zeros(opt.K);
    model.z  = zeros(opt.K, 1);
    model.Z  = zeros(opt.K, opt.N);
    
    % --- Batch processing
    if opt.verbose, before = plotBatchBegin('Init'); end;
    for i=1:ceil(opt.N/batch)
        n1 = (i-1)*batch + 1;
        ne = min(opt.N, i*batch);

        if opt.verbose, before = plotBatch(i, batch, opt.N, 50, before); end;
    
        % Compute subjects log-likelihood
        % -------------------------------
        parfor (n=n1:ne, double(opt.par))
            dat(n) = oneStepLatentInit(dat(n), model, opt);
        end
        
        for n=n1:ne
            
            % Add individual contributions
            % ----------------------------
            model.lb.z.val = model.lb.z.val + dat(n).klz;
            model.Sz       = model.Sz       + dat(n).Sz;
            model.zz       = model.zz       + dat(n).z * dat(n).z';
            model.z        = model.z        + dat(n).z;
            model.Z(:,n)   = dat(n).z;
        end
        
    end
    if opt.verbose, plotBatchEnd; end;

end

% -------------------------------------------------------------------------
%    Batch helper
% -------------------------------------------------------------------------

function before = plotBatchBegin(name)
    before = 0;
    fprintf('%10s | %10s | ', 'Batch', name);
    tic;
end

function before = plotBatch(i, batch, n, total, before)
    exact  = min(n, i*batch)/n*total;
    trunc  = floor(exact);
    step   = trunc - before;
    before = trunc;
    fprintf(repmat('.', 1, step));
end


function plotBatchEnd
    fprintf(' | %fs\n', toc);
end
