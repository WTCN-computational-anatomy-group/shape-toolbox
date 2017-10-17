function [model, dat] = pg_model(opt, dat, model, cont)
% _________________________________________________________________________
%
% FORMAT [model, dat, opt] = pg_registration(opt, (dat), (model))
%
% Optimises a "principal geodesic" model on a collection of images.
%
% The input option structure (opt) should at least contain the field:
% fnames.f - observed images (as a list of filenames)
%
% Alternatively, an input data structure array (dat), of size the number
% of images, can be provided with at least the field:
% f        - observed images (as a list of array or file_array)
%
% The following parameters can be overriden by specifying them in the
% input option (opt) structure:
% model        - Generative model [struct('name', 'normal', 'sigma2', 1)]
% K            - Number of principal geodesics [32]
% lat          - Template lattice dimensions [auto]
% vs           - Template lattice voxel size [auto]
% Mf           - Force same voxel-to-world to all images [read from file]
% itrp         - Interpolation order [1]
% bnd          - Boundary conditions [1]
% emit         - Number of EM iterations [100]
% gnit         - Number of GN iterations for latent variables update [2]
% lsit         - Number of line search iterations [6]
% itgr         - Number of integration steps for geodesic shooting [auto]
% prm          - Parameters of the geodesic differential operator
%                [1e-4 1e-3 0.2 0.05 0.2]
% wpz          - Weights on both parts (A and W'LW) of the prior on z [1 1]
% wpz0         - Initial weights for more robustness [1 5]
% n0           - Number of degrees of freedom of the Wishart prior [20]
% verbose      - Talk during line search [true]
% debug        - Further debuging talk [false]
% loop         - How to split array processing 'subject', 'slice' or 'none'
%                ['subject]
% par          - Parallelise processing 0/n/inf [inf]
% batch        - Batch size for parallelisation [auto].
% directory    - Directory where to store result arrays ['.']
% fnames.result- Filename for the result environment saved after each EM
%                iteration ['pg_result.mat']
% fnames.model - Structure of filenames for all temporary arrays
%                (mu, gmu, (a), w, dw, g, h, ww, A, z, zz, S)
% fnames.dat   - Structure of filenames for all temporary arrays
%                (f, v, ipsi, iphi, pf, c, wmu, gv, hv, z, zz, S)
% ondisk.model - Structure of logical for temporary array
% ondisk.dat   - "      "       "       "       "       "
% _________________________________________________________________________
%
% FORMAT [model, dat, opt] = pg_registration(opt, dat, model, 'continue')
%
% The returned structures (or the saved environment which also contains
% them) can be used as input to start optimising from a previous state.
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
    
    [opt, dat]        = pg_model_input(opt, dat);
    opt               = pg_model_default(opt);
    [opt, dat, model] = pg_model_data(opt, dat, model);
    
    % ---------------------------------------------------------------------
    %    Initialise all variables
    % ---------------------------------------------------------------------
    if ~cont
        if opt.verbose
            fprintf(['%10s | %10s | ' repmat('=',1,50) ' |\n'], 'EM', 'Init');
        end
        [dat, model] = initAll(dat, model, opt);
    end
    
%     model = plotAll(model, opt);
    
    % ---------------------------------------------------------------------
    %    Variable factors
    % ---------------------------------------------------------------------
    % Armijo factor for Subspace line search
    % Because we don't factor part of the log-likelihood when computing the
    % gradient and hessian (in particular, the log(E[ZZ']W'LW) part), it is
    % necessary to check that we don't overshoot. Additionnaly, we try to
    % guess by how much we will overshoot based on the previous EM 
    % iteration.
    model.armijo = 1;
    % In our model ln p(z | W,A) = wpz1 * ln p(z | A) + wpz2 * ln p(z | W)
    % We allow to start with higher or lower weights, and to only use the
    % final weight for the lasts iterations.
    wpzscl1 = logspace(log10(opt.wpz0(1)/opt.wpz(1)), log10(1), opt.emit);
    wpzscl2 = logspace(log10(opt.wpz0(2)/opt.wpz(2)), log10(1), opt.emit);
    
    % ---------------------------------------------------------------------
    %    EM iterations
    % ---------------------------------------------------------------------
    for emit = 1:opt.emit
    
        if opt.verbose
            fprintf(['%10s | %10d | ' repmat('=',1,50) ' |\n'], 'EM', emit);
        end
        
        % Update weights on precision Z
        model.wpz(1) = opt.wpz(1) * wpzscl1(emit);
        model.wpz(2) = opt.wpz(2) * wpzscl2(emit);
        model.regz = model.wpz(1) * model.A + model.wpz(2) * model.ww;
        
        % -----------------------------------------------------------------
        %    M-step (Principal subspace)
        % -----------------------------------------------------------------

        [dat, model] = batchProcess('M', dat, model, opt);
        
        % Factor of the prior : ln p(z|W) + ln p(W)
        % -------------------
        reg = model.wpz(2) * (numeric(model.zz) + numeric(model.S)) ...
              + eye(size(model.zz));
        
        % Gradient
        % --------
        for k=1:opt.K
            lw = spm_diffeo('vel2mom', single(model.w(:,:,:,:,k)), [opt.vs, opt.prm]);
            model.g(:,:,:,:,k) = model.g(:,:,:,:,k) + reg(k,k) * lw;
        end
        
        % Search direction
        % ----------------
        model.dw = prepareOnDisk(model.dw, size(model.w));
        for k=1:opt.K
            model.dw(:,:,:,:,k) = -spm_diffeo('fmg', ...
                single(model.h(:,:,:,:,k)), single(model.g(:,:,:,:,k)), ...
                double([opt.vs reg(k,k) * opt.prm 2 2]));
        end
        model.g = rmarray(model.g);
        model.h   = rmarray(model.h);
        
        [~, model, dat] = lsSubspace(model.dw, model, dat, opt);
        
        model.regz = model.wpz(1) * model.A + model.wpz(2) * model.ww;
        
        % -----------
        % Lower bound
        model.llw = llPriorSubspace(model.w, model.ww, opt.vs, opt.prm);
        model.lba = lbPrecisionZ(model.A, opt.N, model.n0);
        model.lbz = lbLatent(dat, model, opt);
        model     = plotAll(model, opt);
        % -----------
        
        % -----------------------------------------------------------------
        %    E-step (Subject variables)
        % -----------------------------------------------------------------
        
        % Update q(z)
        % -----------
        [dat, model] = batchProcess('E', dat, model, opt);
        
        % -----------
        % Lower bound
        model.llw = llPriorSubspace(model.w, model.ww, opt.vs, opt.prm);
        model.lba = lbPrecisionZ(model.A, opt.N, model.n0);
        model.lbz = lbLatent(dat, model, opt);
        model     = plotAll(model, opt);
        % -----------
        
%         % Centre latent coordinates
%         % -------------------------
%         % > It should accelerate convergence
%         dat = batchProcess('Centre', dat, opt, model.z/opt.N);
%         model.zz = saveOnDisk(model.zz, numeric(model.zz) - model.z * model.z' / opt.N);
%         model.z(:)  = 0;
%         A = precisionZWishart(numeric(model.n0, ...
%             numeric(model.zz) + numeric(model.S), opt.N);
%         model.A = saveOnDisk(model.A, A);
%         
%         % -----------
%         % Lower bound
%         model.llw = llPriorSubspace(model.w, model.ww, opt.vs, opt.prm);
%         model.lba = lbPrecisionZ(model.A, opt.N, model.n0);
%         model.lbz = lbLatent(dat, model, opt);
%         model     = plotAll(model, opt);
%         % -----------
        
        % Orthogonalise
        % -------------
        if opt.verbose, fprintf('%10s | %10s ', 'Ortho', ''); tic; end;
        % > Needed for our Hessian approximation 
        [U, iU] = orthogonalisationMatrix(numeric(model.zz), numeric(model.ww));
        if opt.verbose, fprintf('| %6gs\n', toc); end;
        
        % Rescale
        % -------
        if opt.verbose, fprintf('%10s | %10s ', 'Rescale', ''); tic; end;
        ezz = U*(numeric(model.zz)+numeric(model.S))*U';
        if model.n0 == 0
            [Q, iQ] = scalePG(opt.N, opt.K);
        else
            [Q, iQ] = gnScalePG(ezz, model.n0, opt.N, model.wpz(2));
        end
        if opt.verbose, fprintf('| %6gs\n', toc); end;
        Q = Q*U;
        iQ = iU*iQ;
        [model, dat] = rotateAll(model, dat, opt, Q, iQ);
        A = precisionZWishart(model.n0, ...
            numeric(model.zz) + numeric(model.S), opt.N);
        model.A = saveOnDisk(model.A, A);
        
        model.regz = model.wpz(1) * model.A + model.wpz(2) * model.ww;
        
        % -----------
        % Lower bound
        model.llw = llPriorSubspace(model.w, model.ww, opt.vs, opt.prm);
        model.lba = lbPrecisionZ(model.A, opt.N, model.n0);
        model.lbz = lbLatent(dat, model, opt);
        model     = plotAll(model, opt);
        % -----------
        
        % -----------------------------------------------------------------
        %    Update template
        % -----------------------------------------------------------------
        if opt.verbose, fprintf('%10s | %10s ', 'Template', ''); tic; end;
        if opt.tpm
            model.a = updateMuML(opt.model, dat.pf, dat.c, ...
                                 'par', opt.par, 'debug', opt.debug, ...
                                 'output', model.a);
            model.gmu = templateGrad(model.a, opt.itrp, opt.bnd, ...
                'debug', opt.debug, 'output', model.gmu);
            model.mu = reconstructProbaTemplate(model.a, ...
                'loop', '', 'par', opt.par, 'debug', opt.debug, ...
                'output', model.mu);
        else
            model.mu = updateMuML(opt.model, dat.pf, dat.c, ...
                                  'par', opt.par, 'debug', opt.debug, ...
                                  'output', model.mu);
            model.gmu = templateGrad(model.mu, opt.itrp, opt.bnd, ...
                'debug', opt.debug, 'output', model.gmu);
        end
        if opt.verbose, fprintf('| %6gs\n', toc); end;
        
        % -----------
        % Lower bound
        dat = batchProcess('Update', dat, model, opt, 'llm');
        model.llm = 0;
        for n=1:opt.N
            model.llm = model.llm + dat(n).llm;
        end
        model     = plotAll(model, opt);
        % -----------
        
        save(fullfile(opt.directory, opt.fnames.result), 'model', 'dat', 'opt');
        
    end % < EM loop
    
end

% =========================================================================

function model = plotAll(model, opt)
    if opt.verbose
        % Template & PG
        subplot(3, 3, 1)
        imagesc(model.mu(:,:,ceil(size(model.mu,3)/2),1));
        title('template')
%         subplot(3, 3, 2)
%         imagesc(model.w(:,:,ceil(size(model.mu,3)/2),1,1));
%         title('PG1 x')
        subplot(3, 3, 2)
        imagesc(model.w(:,:,ceil(size(model.mu,3)/2),2,1));
        title('PG1 y')
        % Precision
        subplot(3, 3, 4)
        imagesc(model.ww)
        colorbar
        title('E*[W''LW]')
        subplot(3, 3, 5)
        imagesc(model.A)
        colorbar
        title('E*[A]')
        % Lower bound
        subplot(3, 3, 3)
        if ~isfield(model, 'lb')
            model.lb = [];
        end
        model.lb = [model.lb (model.llm + model.llw + model.lba + model.lbz)];
        plot(model.lb)
        title('Lower bound')
        subplot(3, 3, 6)
        if ~isfield(model, 'savellm')
            model.savellm = [];
        end
        model.savellm = [model.savellm model.llm];
        plot(model.savellm, 'r-')
        title('Matching part')
        subplot(3, 3, 7)
        if ~isfield(model, 'savellw')
            model.savellw = [];
        end
        model.savellw = [model.savellw model.llw];
        plot(model.savellw, 'g-')
        title('Subspace part')
        subplot(3, 3, 8)
        if ~isfield(model, 'savelbz')
            model.savelbz = [];
        end
        model.savelbz = [model.savelbz model.lbz];
        plot(model.savelbz, 'k-')
        title('Latent part')
        subplot(3, 3, 9)
        if ~isfield(model, 'savelba')
            model.savelba = [];
        end
        model.savelba = [model.savelba model.lba];
        plot(model.savelba, 'c-')
        title('Precision part')
        
        drawnow
        
        fprintf('%10s | ', 'LB');
        if length(model.lb) > 1
            if model.lb(end) - model.lb(end-1) > 0
                fprintf('%10s | ', '(+)');
            elseif model.lb(end) - model.lb(end-1) < 0
                fprintf('%10s | ', '(-)');
            else
                fprintf('%10s | ', '(=)');
            end
        else
            fprintf('%10s | ', '');
        end
        fprintf(' %6g\n', model.lb(end));
    end
end

function [dat, model] = initAll(dat, model, opt)

    % Init identity transforms
    % ------------------------
    
    % --- Zero init of W
    if ~checkarray(model.w)
        model.w = initSubspace(opt.lat, opt.K, 'type', 'zero', ...
            'debug', opt.debug, 'output', model.w);
    end
    if ~checkarray(model.ww)
        model.ww = eye(opt.K);
%         model.ww = precisionZ(model.w, opt.vs, opt.prm, ...
%             'debug', opt.debug, 'output', model.ww) + eye(opt.K);
    end
    
    % --- Zero init of Z
    [dat, model] = batchProcess('Init', 'zero', dat, model, opt);
    
    % --- Init of subject specific arrays
    dat = batchProcess('Update', dat, model, opt, ...
        {'v', 'ipsi', 'iphi', 'pf', 'c'}, 'clean', {'ipsi', 'iphi'});

    % --- Init template + Compute template spatial gradients + Build TPMs
    if opt.tpm
        model.a = updateMuML(opt.model, dat.pf, dat.c, ...
                             'par', opt.par, 'debug', opt.debug, ...
                             'output', model.a);
        model.gmu = templateGrad(model.a, opt.itrp, opt.bnd, ...
            'debug', opt.debug, 'output', model.gmu);
        model.mu = reconstructProbaTemplate(model.a, ...
            'loop', '', 'par', opt.par, 'debug', opt.debug, ...
            'output', model.mu);
    else
        model.mu = updateMuML(opt.model, dat.pf, dat.c, ...
                              'par', opt.par, 'debug', opt.debug, ...
                              'output', model.mu);
        model.gmu = templateGrad(model.mu, opt.itrp, opt.bnd, ...
            'debug', opt.debug, 'output', model.gmu);
    end
    
    % Compute initial matching LL
    % ---------------------------
    dat = batchProcess('Update', dat, model, opt, {'llm'});
    model.llm = 0;
    for n=1:opt.N
        model.llm = model.llm + dat(n).llm;
    end
    
    % Init latent coordinates
    % -----------------------

    % --- Random init of E[z]
    [dat, model] = batchProcess('Init', 'rand', dat, model, opt);

    % --- Orthogonalise sum{E[z]E[z]'}
    [U,S] = svd(model.zz);
    Rz    = 0.1*sqrt(opt.N/opt.K)*U/diag(sqrt(diag(S)+eps));
    dat   = batchProcess('rotate', dat, opt, Rz');
    model.zz = Rz' * model.zz * Rz;
    model.S  = Rz' * model.S * Rz;
    
    % Init precision of z
    % -------------------
    model.A = saveOnDisk(model.A, eye(opt.K));

end
