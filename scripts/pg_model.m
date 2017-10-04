function [model, dat] = pg_model(opt, dat, model, cont)
% _________________________________________________________________________
%
% FORMAT [model, dat, opt] = pg_registration(opt, (dat), (model))
%
% Optimizes a "principal geodesic" model on a collection of images.
%
% The input option structure (opt) should at least contain the field:
% fnames.f - observed images (as a list of filenames)
%
% Alternatively, an input data structure array (dat), of size the number
% of images, can be provided with at least the field:
% f        - observed images (as an array or file_array)
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
%                (mu, gmu, (a), w, dw, g, h, ww, A0, A, z, zz, S)
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
        [dat, model] = initAll(dat, model, opt);
    end
    
    plotAll(model, opt)
    
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
    % In our model ln p(z | W,A) = ln p(z | A) + wpz * ln p(z | W)
    % ln p(z | W) is necessary during the first iterations, when the
    % Wishart prior is not very accurate, as a "smoothing" constraint.
    % Its weight is lowered after each EM iteration.
    % The rescaling starts with 0.5 and the cumulative product is very
    % quickly close to 0.
    wpzscl = logspace(log10(0.5), 0, opt.emit);
    
    % ---------------------------------------------------------------------
    %    EM iterations
    % ---------------------------------------------------------------------
    for emit = 1:opt.emit
    
        if opt.verbose
            fprintf(['EM %4d | ' repmat('=',1,70) '\n'], emit);
        end
        
        % -----------------------------------------------------------------
        %    M-step (Principal subspace)
        % -----------------------------------------------------------------

        [dat, model] = batchProcess('M', dat, model, opt);
        
        % Factor of the prior : ln p(z|W) + ln p(W)
        % -------------------
        reg = model.wpz(2) * numeric(model.zz) + model.wpz(1) * eye(size(model.zz));
        
        % Gradient
        % --------
        for k=1:opt.K
            lw = spm_diffeo('vel2mom', single(model.w(:,:,:,:,k)), [opt.vs, opt.prm]);
            model.g(:,:,:,:,k) = model.g(:,:,:,:,k) + reg(k,k) * lw;
        end
        
        % Search direction
        % ----------------
        for k=1:opt.K
            model.dw(:,:,:,:,k) = -spm_diffeo('fmg', ...
                single(model.h(:,:,:,:,k)), single(model.g(:,:,:,:,k)), ...
                double([opt.vs reg(k,k) * opt.prm 2 2]));
        end
        model.g = rmarray(model.g);
        model.h = rmarray(model.h);
        
        [~, model, dat] = lsSubspace(model.dw, model, dat, opt);
        
        plotAll(model, opt);
        
        % -----------------------------------------------------------------
        %    E-step (Subject variables)
        % -----------------------------------------------------------------
        
        [dat, model] = batchProcess('E', dat, model, opt);
        
        % -----------------------------------------------------------------
        %    Update template
        % -----------------------------------------------------------------
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
        
        plotAll(model, opt);
        
        save(fullfile(opt.directory, opt.fnames.result), 'model', 'dat', 'opt');
        
        model.wpz(2) = model.wpz(2) * wpzscl(emit);
        
    end % < EM loop
    
end

% =========================================================================

function plotAll(model, opt)
    if opt.verbose
        subplot(2, 3, 1)
        imagesc(model.mu(:,:,ceil(size(model.mu,3)/2),1));
        title('template')
        subplot(2, 3, 2)
        imagesc(model.w(:,:,ceil(size(model.mu,3)/2),1,1));
        title('PG1 x')
        subplot(2, 3, 3)
        imagesc(model.w(:,:,ceil(size(model.mu,3)/2),2,1));
        title('PG1 y')
        subplot(2, 3, 4)
        imagesc(model.ww)
        colorbar
        title('W''LW')
        subplot(2, 3, 5)
        imagesc(model.zz+model.S)
        colorbar
        title('E[zz'']')
        subplot(2, 3, 6)
        imagesc(model.A)
        colorbar
        title('A')
        drawnow
    end
end

function [dat, model] = initAll(dat, model, opt)

    % Init latent coordinates
    % -----------------------

    % --- Zero init of W
    if ~checkarray(model.w)
        model.w = initSubspace(opt.lat, opt.K, 'type', 'zero', ...
            'debug', opt.debug, 'output', model.w);
    end
    if ~checkarray(model.ww)
        model.ww = precisionZ(model.w, opt.vs, opt.prm, ...
            'debug', opt.debug, 'output', model.ww);
    end

    % --- Random init of E[z]
    [dat, model] = batchProcess('init', dat, model, opt);

    % --- Orthogonalise sum{E[z]E[z]'}
    [U,S] = svd(model.zz);
    Rz    = 0.1*sqrt(opt.N/opt.K)*U/diag(sqrt(diag(S)+eps));
    dat   = batchProcess('rotate', dat, opt, Rz);
    model.zz = Rz' * model.zz * Rz;
    model.S  = Rz' * model.S * Rz;

    % Init precision of z
    % -------------------
    model.A = saveOnDisk(model.A, numeric(model.A0) * model.n0);
    model.n = model.n0;

    % Init subject-specific arrays
    % ----------------------------
    dat = batchProcess('Update', dat, model, opt, ...
        {'v', 'ipsi', 'iphi', 'pf', 'c'}, 'clean', {'ipsi', 'iphi'});

    % Init template + Compute template spatial gradients + Build TPMs
    % ---------------------------------------------------------------
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

    % Start with M step
    % -----------------
    % It is easier to generate random latent coordinate than random
    % principal components.

    [dat, model] = stepM0(dat, model, opt);
    plotAll(model, opt)
    [dat, model] = batchProcess('E', dat, model, opt);

    %  Update template
    % ----------------
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

end

function [dat, model] = stepM0(dat, model, opt)

    % Init log-likelihood
    % -------------------
    model.llz  = 0;
    model.lllz = 0;
    model.llw  = 0;
    model.ww   = saveOnDisk(model.ww, eye(opt.K));
    
    % Matchign term
    % -----------------
    [dat, model] = batchProcess('M0', dat, model, opt);

    % Factor of the prior : ln p(z|W) + ln p(W)
    % -------------------
    reg = model.wpz(2) * numeric(model.zz) + model.wpz(1) * eye(size(model.zz));
        
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
    
    [~, model, dat] = lsSubspace(model.dw, model, dat, opt);
end