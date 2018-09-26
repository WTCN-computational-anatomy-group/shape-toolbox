function [model, dat, opt] = pgva_model(varargin)
% _________________________________________________________________________
%
%     PGVA model (2018) - Principal Geodesic from Velocity + Affine
% _________________________________________________________________________
%
% FORMAT [model, dat, opt] = pgva_model(input, (opt))
%
% Learn a principal subspace of deformation from data and/or velocity 
% fields.
%
% This model relies on modelling velocity fields (which encode non-linear
% diffeomorphic transforms) as
%       v = W * z + r.
% In our structure, we will consequently name
%   v    = velocity field
%   w    = principal subspace of deformation
%   z    = latent coordinates in the principal subspace
%   r    = residual field
%   q    = parameters of a rigid (or affine) transform, to align shapes
%   a/mu = Template (i.e., mean shape)
%
% In this model, we alternate between explicitely fitting velocity fields
% (v) by Gauss-Newton optimisation and computing closed-form solutions of 
% the principal geodesic decomposition (W and z).
%
% -------------------------------------------------------------------------
%
% MANDATORY INPUT FILES
% ---------------------
% The input option structure (opt) should contain at least one of the 
% fields:
% input.f - observed data images     (as a list of filenames)
% input.v - observed velocity fields (as a list of filenames)
%
% OPTIONAL INPUT FILES
% --------------------
% Additionnaly, starting estimates for some model parameters can be
% provided:
% input.w  - principal subspace
% input.a  - log-template (for categorical cases)
% input.mu - template (for intensity cases)
%
% NB
% --
% The principal subspace, template and velocity fields should all have
% identical dimensions and voxel sizes.
% Input images can have different dimensions.
% 
% -------------------------------------------------------------------------
% The following parameters can be overriden by specifying them in the
% input option (opt) structure:
%
% MODEL
% -----
% model.name   - Generative data model ['normal']/'categorical'/'bernoulli'
% model.sigma2 - If normal model: initial noise variance estimate [1]
% model.nc     - (categorical only) Number of classes [from input]
% pg.K     - Number of principal geodesics [32]
% pg.prm   - Parameters of the geodesic operator [0.001 0 10 0.1 0.2]
% pg.bnd   - Boundary conditions for the geodesic operator [1 = circulant]
% mixreg.a0- Prior expected value of the mixture weight [0.5]
% mixreg.n0- Prior DF of the mixture weight [1e-4]
% tpl.vs   - Lattice voxel size [auto]
% tpl.lat  - Lattice dimensions [auto]
% tpl.prm  - Parameters of the field operator [1e-3  1e-1 0]
% tpl.bnd  - Boundary conditions for the field operator [1 = circulant]
% tpl.itrp - Interpolation order [1]
% tpl.ld   - Field operator log-determinant [from prm]
% v.l0     - Prior expected anatomical noise precision [17]
% v.n0     - Prior DF of the anatomical noise precision [10]
% z.init   - Latent initialisation mode ['auto']/'zero'/'rand'
% z.A0     - Prior expected latent precision matrix [eye(K)] 
% z.n0     - Prior DF of the latent precision matrix [K]
% q.A0     - Prior expected affine precision matrix [eye(M)]
% q.n0     - Prior DF of the affine precision matrix [M]
% q.B      - Affine_basis ['rigid']
% q.hapx   - Approximate affine hessian [true]
% f.M      - Force same voxel-to-world to all images [read from file]
%
% optimise.pg.w      - Optimise subspace                [true]
% optimise.z.z       - Optimise latent coordinates      [true]
% optimise.z.A       - Optimise latent precision        [true]
% optimise.q.q       - Optimise affine coordinates      [true]
% optimise.q.A       - Optimise affine precision        [true]
% optimise.v.v       - Optimise velocity fields         [true]
% optimise.v.l       - Optimise residual precision      [true]
% optimise.tpl.a     - Optimise template                [true]
% optimise.tpl.scale - Optimise template scaling        [false]
% optimise.mixreg.w  - Optimise mixture weight          [true]
% optimise.mixreg.a  - Optimise mixture weight prior    [true]
%
% PROCESSING
% ----------
% match        - Matching term version 'push'/['pull']
% iter.em      - Maximum number of EM iterations [1000]
% iter.gn      - Maximum number of Gauss-Newton iterations [1]
% iter.ls      - Maximum number of line search iterations [6]
% iter.itg     - Number of integration steps for geodesic shooting [auto]
% iter.pena    - Penalise Gauss-Newton failures [true]
% lb.threshold - Convergence criterion (lower bound gain) [1e-5]
% lb.moving    - Moving average over LB gain [3]
% split.loop   - How to split array processing: 'none'/'slice'/['subject']
% split.par    - Parallelise processing (number of workers): 0/n/[inf]
% split.batch  - Batch size for parallelisation [auto]
% ui.verbose   - Talk during processing [true]
% ui.debug     - Further debuging talk [false]
% ui.ftrack    - Figure object for the lower bound tracking [gcf]
% dist         - Distributed processing. See `help distribute_default`.
%
% I/O
% ---
% dir.model     - Directory where to store model arrays and workspace ['.']
% dir.dat       - Directory where to store data array [next to input]
% fnames.result - Filename for the result environment saved after each EM
%                 iteration ['pgva_model.mat']
% fnames.model  - Structure of filenames for all file arrays
% fnames.dat    - Structure of filenames for all file arrays
% ondisk.model  - Structure of logical for temporary array [default_ondisk]
% ondisk.dat    - "      "       "       "       "       "
% _________________________________________________________________________
%
% FORMAT [model, dat, opt] = pgva_model(opt, dat)
%
% The returned structures (or the saved environment which also contains
% them) can be used as input to start optimising from a previous state.
% _________________________________________________________________________

% _________________________________________________________________________
%
%                             Graphical model
% =========================================================================
% 
% Velocity part
% -------------
%                          [L] --> (w)
%                                   |
%                                   v
%  [Az0, nz0] --> (Az) --> (z) --> (v) <-- (lam) <-- [lam0, nl0]
%                                   ^
%                                   |
%                                  [L]
%
% Data part
% ---------
%
% (v) --> <phi> --> <psi> <-- <xi> <-- <R> <-- (q) <-- (Aq) <-- [Aq0, nq0]
%                     |
%                     v
%   [La] --> (a) --> <mu> --> {f}
%
% Legend
% ------
% { } = observed variable
% ( ) = latent variable
% < > = deterministic variable
% [ ] = fixed parameter
% _________________________________________________________________________

    cleanupObj = hello('PGVA model');
    
    % -----------
    % Parse input
    % -----------
    if nargin == 0
        help pgva_model
        error('At least one input argument is needed.')
    end
    if nargin >= 3
        opt   = varargin{1};
        dat   = varargin{2};
        model = varargin{3};
        cont  = true;
    else
        input = varargin{1};
        if nargin >= 2
            opt = varargin{2};
        else
            opt = struct;
        end
        cont = false;
    end
    
    % -----------
    % Add to path
    % -----------
    setpath('pgva');
    
    % ---------------------------------------------------------------------
    %    Initialisation
    % ---------------------------------------------------------------------
    if ~cont
        
        % Default parameters + prepare structures & file arrays
        % -----------------------------------------------------------------
        [opt,dat,model] = pgva_model_input(input,opt);    % Read observed
        opt             = pgva_model_default(opt);        % Read options
        [opt,dat,model] = pgva_model_data(opt,dat,model); % Set arrays

        % Copy screen output to file
        % --------------------------
        if ~isempty(opt.fnames.log)
            if exist(fullfile(opt.dir.model, opt.fnames.log), 'file')
                delete(fullfile(opt.dir.model, opt.fnames.log));
            end
            diary(fullfile(opt.dir.model, opt.fnames.log));
        end
        
        % Post-set parameters
        % -----------------------------------------------------------------
        % We store some values to avoid unneeded computation
        spm_diffeo('boundary', opt.pg.bnd);
%         [~, opt.pg.ld] = spm_shoot_greens('kernel', double(opt.tpl.lat), double([opt.tpl.vs opt.pg.prm]), opt.pg.bnd);
%         opt.pg.ld = opt.pg.ld(1);
        opt.pg.ld = ldapprox(opt.pg.prm, 'vs', opt.tpl.vs, 'dim', [opt.tpl.lat 3]);
        ker = spm_diffeo('kernel', double(opt.tpl.lat), double([opt.tpl.vs opt.pg.prm]));
        opt.pg.ker = [ker(1,1,1,1,1) ker(1,1,1,2,2) ker(1,1,1,3,3)];
        clear ker
        
        % Initialise all arrays (= model variables)
        % -----------------------------------------------------------------
        if opt.ui.verbose
            fprintf(['%10s | %10s | ' repmat('=',1,50) ' |\n'], 'EM', 'Init');
        end
        [dat, model] = pgva_model_init(dat, model, opt);
        
        % Some more stuff regarding processing
        % -----------------------------------------------------------------
        model.emit      = 1;
        model.q.active  = true;
        model.v.active  = true;
        model.pg.active = true;
    end
    
    % ---------------------------------------------------------------------
    %    EM iterations
    % ---------------------------------------------------------------------
    for emit = model.emit:opt.iter.em
        
        % -----------------------------------------------------------------
        %    Lower bound gain
        % -----------------------------------------------------------------
        model = updateLowerBound(model, 'gain');
        
        % -----------------------------------------------------------------
        %    Save current state
        % -----------------------------------------------------------------
        if ~isempty(opt.fnames.result)
            % Ensure nifti headers are ok
            createAllNifti(dat, model, opt);
            % Write workspace
            ftrack = opt.ui.ftrack;
            opt.ui.ftrack = nan;
            save(fullfile(opt.dir.model, opt.fnames.result), ...
                 'model', 'dat', 'opt');
            opt.ui.ftrack = ftrack;
        end
        pgva_plot_all(model, opt);
        if ~isempty(opt.fnames.fig) && (~islogical(opt.ui.ftrack) || opt.ui.ftrack)
            saveas(gcf, fullfile(opt.dir.model, opt.fnames.fig));
        end
        
        % -----------------------------------------------------------------
        %    General tracking
        % -----------------------------------------------------------------
        
        if opt.ui.verbose
            fprintf(['%10s | %10d | ' repmat('=',1,50) ' |\n'], 'EM', emit);
        end
        
        N = numel(model.lb.lb.gainlist);
        moving_gain = mean(abs(model.lb.lb.gainlist(N:-1:max(1,N-opt.lb.moving+1))));
        if moving_gain < opt.lb.threshold
            if opt.optimise.q.q && ~model.q.active
                model.q.active = true;
                fprintf('%10s | %10s\n', 'Activate', 'Affine');
            elseif opt.optimise.v.v && ~model.v.active
                model.v.active = true;
                fprintf('%10s | %10s\n', 'Activate', 'Velocity');
            elseif (opt.optimise.pg.w || opt.optimise.z.z) && ~model.pg.active
                model.pg.active = true;
                fprintf('%10s | %10s\n', 'Activate', 'PG');
            else
                fprintf('Converged :D\n');
                break
            end
        end
        model.emit = emit;
        
        % -----------------------------------------------------------------
        %    Affine
        % -----------------------------------------------------------------
        if model.q.active && opt.f.N
        
            if opt.optimise.q.q
                % Update parameters
                % -----------------
                [dat, model] = pgva_batch('FitAffine', dat, model, opt);

                % -----------
                % Lower bound
                model           = updateLowerBound(model);
                pgva_plot_all(model, opt);
                % -----------
            end
            
            if opt.optimise.q.A && opt.q.Mr
                % Update prior
                % ------------
                rind = opt.q.rind;
                model.q.A = spm_prob('Wishart', 'up', ...
                                     opt.f.N, 0, model.q.qq(rind,rind) + model.q.S(rind,rind), ...
                                     opt.q.A0, opt.q.n0);

                % -----------
                % Lower bound
                [dat, model]    = pgva_batch('LB', 'PrecisionQ', dat, model, opt);
                if opt.q.n0
                    model.lb.Aq.val = -spm_prob('Wishart', 'kl', ...
                                           model.Aq,         opt.nq0+opt.f.N, ...
                                           eye(numel(rind)), opt.nq0, ...
                                           'normal');
                end
                model = updateLowerBound(model);
                pgva_plot_all(model, opt);
                % -----------
            end
            
        end
        
        % -----------------------------------------------------------------
        %    Velocity field
        % -----------------------------------------------------------------
        if model.v.active
            
            % Update subjects
            % ---------------
            if opt.optimise.v.v && opt.f.N
                [dat, model] = pgva_batch('FitVelocity', dat, model, opt);

                % -----------
                % Lower bound
                model = updateLowerBound(model);
                pgva_plot_all(model, opt);
                % -----------
            end
            
            % Update precision
            % ----------------
            if opt.optimise.v.l
                old_l = model.v.l;
                K = 3*prod(opt.tpl.lat);
                model.v.n = opt.v.n0 + model.mixreg.w(1) * (opt.v.N+opt.f.N);
                model.v.l = opt.v.n0/opt.v.l0 + (model.mixreg.w(1)/K)*(model.v.uncty + model.v.tr + model.v.reg);
                model.v.l = model.v.n/model.v.l;
                if opt.ui.verbose, fprintf('%10s | %10s | %8.6g -> %8.6g\n', 'Lambda', '', old_l, model.v.l); end

                % -----------
                % Lower bound
                [dat, model]    = pgva_batch('LB', 'Lambda', dat, model, opt);
                model = updateLowerBound(model);
                pgva_plot_all(model, opt);
                % -----------
            end
            
        end
        
        % -----------------------------------------------------------------
        %    Principal Geodesic Analysis
        % -----------------------------------------------------------------
        if model.pg.active

            % Principal subspace
            % -------------------------------------------------------------
            if opt.optimise.pg.w
                
                if opt.ui.verbose, fprintf('%10s | %10s ', 'PG', ''); tic; end
                M = model.z.S + model.z.zz + model.pg.n * eye(opt.pg.K) / (model.v.l*model.mixreg.w(1));
                M = spm_matcomp('Inv', M);
                P = model.z.Z' * M;
                for k=1:opt.pg.K
                    w1 = zeros([opt.tpl.lat 3], 'single');
                    for n=1:numel(dat)
                        w1 = w1 + numeric(dat(n).v.v) * P(n,k);
                    end
                    model.pg.w(:,:,:,:,k) = w1;
                end
                clear M P
                model.pg.ww = precisionZ(model.pg.w, opt.tpl.vs, opt.pg.prm);
                if opt.ui.verbose, fprintf('| %6gs\n', toc); end

                % -----------
                % Lower bound
                [dat, model] = pgva_batch('LB', 'Subspace', dat, model, opt);
                model        = updateLowerBound(model);
                pgva_plot_all(model, opt);
                % -----------
                
            end

            % Latent coordinates
            % -------------------------------------------------------------
            if opt.optimise.z.z
                % Coordinates (variational update)
                % --------------------------------
                [dat, model] = pgva_batch('FitLatent', dat, model, opt);

                % -----------
                % Lower bound
                [dat, model] = pgva_batch('LB', 'Latent', dat, model, opt);
                model        = updateLowerBound(model);
                pgva_plot_all(model, opt);
                % -----------
            end

            % Orthogonalisation
            % -------------------------------------------------------------
            [~,p] = chol(model.pg.ww); % Check pos-def
            if opt.optimise.z.z && opt.optimise.pg.w && opt.optimise.z.A && p == 0
                
                % Orthogonalise
                % -------------
                if opt.ui.verbose, fprintf('%10s | %10s ', 'Ortho', ''); tic; end
                [U, iU] = orthogonalisationMatrix(model.z.zz+model.z.S, model.pg.ww);
                if opt.ui.verbose, fprintf('| %6.3fs\n', toc); end

                % Rescale
                % -------
                if opt.ui.verbose, fprintf('%10s | %10s ', 'Rescale', ''); tic; end
                [Q, iQ] = gnScalePG(iU' * model.pg.ww * iU * model.pg.n, ...
                                    U   * model.z.zz  * U', ...
                                    U   * model.z.S   * U', ...
                                    opt.z.A0, opt.z.n0, model.pg.n);
                if opt.ui.verbose, fprintf('| %6.3fs\n', toc); end
                Q  = Q  * U;
                iQ = iU * iQ;

                % Rotate
                % ------
                [dat, model] = pgva_batch('RotateSubspace', Q, iQ, dat, model, opt);
           
                % Latent precision
                % ----------------
                model.z.A = spm_prob('Wishart', 'up', ...
                                     opt.v.N+opt.f.N, 0, model.z.zz + model.z.S, ...
                                     opt.z.A0, opt.z.n0);

                % -----------
                % Lower bound
                [dat, model] = pgva_batch('LB', 'Orthogonalise', dat, model, opt);
                model        = updateLowerBound(model);
                pgva_plot_all(model, opt);
                % -----------
                

            % Latent precision
            % -------------------------------------------------------------
            elseif opt.optimise.z.A
                
                model.z.A = spm_prob('Wishart', 'up', ...
                                     opt.v.N+opt.f.N, 0, model.z.zz + model.z.S, ...
                                     opt.z.A0, opt.z.n0);

                % -----------
                % Lower bound
                [dat, model] = pgva_batch('LB', 'PrecisionZ', dat, model, opt);
                model        = updateLowerBound(model);
                pgva_plot_all(model, opt);
                % -----------
            end
            
        end
        
        % -----------------------------------------------------------------
        %    Template
        % -----------------------------------------------------------------
        if opt.f.N && opt.optimise.tpl.a
            if opt.ui.verbose, fprintf('%10s | %10s ', 'Template', ''); tic; end
            if opt.tpl.cat
                model.tpl.a = updateMuML(opt.model, dat, ...
                                         'lat',    opt.tpl.lat,   ...
                                         'par',    opt.split.par, ...
                                         'debug',  opt.ui.debug,  ...
                                         'output', model.tpl.a);
                model.tpl.gmu = templateGrad(model.tpl.a,  ...
                                             opt.tpl.bnd,  ...
                                             'debug',  opt.ui.debug, ...
                                             'output', model.tpl.gmu);
                model.tpl.mu = reconstructProbaTemplate(model.tpl.a, ...
                                                        'par',    opt.split.par, ...
                                                        'debug',  opt.ui.debug,  ...
                                                        'output', model.tpl.mu);
            else
                model.tpl.mu = updateMuML(opt.model, dat, ...
                                          'lat',    opt.tpl.lat,   ...
                                          'par',    opt.split.par, ...
                                          'debug',  opt.ui.debug,  ...
                                          'output', model.tpl.mu);
                model.tpl.gmu = templateGrad(model.tpl.mu, ...
                                             opt.tpl.bnd, ...
                                             'debug',  opt.ui.debug, ...
                                             'output', model.tpl.gmu);
            end
            if opt.ui.verbose, fprintf('| %6.3s\n', toc); end
            % -----------
            % Lower bound
            [dat, model] = pgva_batch('LB', 'Template', dat, model, opt);
            model        = updateLowerBound(model);
            pgva_plot_all(model, opt);
            % -----------
        end
        
        % -----------------------------------------------------------------
        %    Mixture
        % -----------------------------------------------------------------
        if opt.optimise.mixreg.w || opt.optimise.mixreg.a
            if opt.optimise.mixreg.w
                old_w = model.mixreg.w(1);
                N = model.mixreg.n;
                A = model.mixreg.a;
                dg = spm_prob('DiGamma', N*A) - spm_prob('DiGamma', N);
                rho1 = dg;
                rho2 = -dg;
                for n=1:numel(dat)
                    rho1 = rho1 + dat(n).v.lb.ll1;
                    rho2 = rho2 + dat(n).v.lb.ll2;
                end
                model.mixreg.w(1) = 1/(1+exp(rho2-rho1));
                model.mixreg.w(2) = 1 - model.mixreg.w(1);
                if opt.ui.verbose, fprintf('%10s | %10s | %8.6g -> %8.6g\n', 'Mix W', '', old_w, model.mixreg.w(1)); end
            end
            if opt.optimise.mixreg.a
                old_a = model.mixreg.a;
                N0 = opt.mixreg.n0;
                A0 = opt.mixreg.a0;
                model.mixreg.a = (N0*A0 + model.mixreg.w(1))/(N0+1);
                if opt.ui.verbose, fprintf('%10s | %10s | %8.6g -> %8.6g\n', 'Mix A', '', old_a, model.mixreg.a); end
            end
            N = model.mixreg.n;
            A = model.mixreg.a;
            dg = spm_prob('DiGamma', N*A) - spm_prob('DiGamma', N);
            dig = spm_prob('DiGamma', N*(1-A)) - spm_prob('DiGamma', N);
            if opt.optimise.mixreg.w
                W1 = model.mixreg.w(1);
                W2 = model.mixreg.w(2);
                model.lb.wr.val = W1*(log(W1+eps) - dg) + W2*(log(W2+eps) + dg);
                model.lb.wr.type = '-kl';
            end
            if opt.optimise.mixreg.a
                model.lb.ar.val =   (N*A - N0*A0) * dg ...
                                  + (N*(1-A) - N0*(1-A0)) * dig...
                                  + beta_norm(N0*A0,N0*(1-A0)) ...
                                  - beta_norm(N*A,N*(1-A));
                model.lb.ar.type = '-kl';
            end
            % -----------
            % Lower bound
            [dat, model] = pgva_batch('LB', 'Mixture', dat, model, opt);
            model = updateLowerBound(model);
            pgva_plot_all(model, opt);
            % -----------
        end
    end
end

% =========================================================================
function b = beta_norm(a,b)
    b = gammaln(a) + gammaln(b) - gammaln(a+b);
end