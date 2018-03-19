function [model, dat, opt] = pgra_model(varargin)
% _________________________________________________________________________
%
%     PGRA model (2017) - Principal Geodesic + Residual + Affine
% _________________________________________________________________________
%
% FORMAT [model, dat, opt] = pgra_model(input, (opt))
%
% Learn a principal subspace of deformation from data.
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
% In this model, W, z, and r are explicit in the model and are all  
% obtained by Gauss-Newton optimisation.
%
% -------------------------------------------------------------------------
%
% MANDATORY INPUT FILES
% ---------------------
% The input option structure (opt) should contain at least one of the 
% fields:
% input.f - observed data images     (as a list of filenames)
%
% OPTIONAL INPUT FILES
% --------------------
% Additionnaly, starting estimates for some model parameters can be
% provided:
% input.w  - principal subspace
% input.a  - log-template (for categorical cases: categorical or bernoulli)
% input.mu - template     (for intensity cases: normal or laplace)
%
% NB
% --
% The principal subspace and template should have identical dimensions and 
% voxel sizes.
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
% pg.K     - Number of principal geodesics [32]
% pg.prm   - Parameters of the geodesic operator [1e-4 1e-3 0.2 0.05 0.2]
% pg.bnd   - Boundary conditions for the geodesic operator [1 = circulant]
% mixreg.a0- Prior expected value of the mixture weight [0.5]
% mixreg.n0- Prior DF of the mixture weight [1e-4]
% tpl.vs   - Lattice voxel size [auto]
% tpl.lat  - Lattice dimensions [auto]
% tpl.prm  - Parameters of the field operator [1e-3  1e-1 0]
% tpl.bnd  - Boundary conditions for the field operator [1 = circulant]
% tpl.itrp - Interpolation order [1]
% tpl.ld   - Field operator log-determinant [from prm]
% r.l0     - Prior expected anatomical noise precision [17]
% r.n0     - Prior DF of the anatomical noise precision [10]
% z.init   - Latent initialisation mode ['auto']/'zero'/'rand'
% z.A0     - Prior expected latent precision matrix [eye(K)] 
% z.n0     - Prior DF of the latent precision matrix [K]
% q.A0     - Prior expected affine precision matrix [eye(M)]
% q.n0     - Prior DF of the affine precision matrix [M]
% q.B      - Affine_basis ['rigid']
% q.hapx   - Approximate affine hessian [true]
% f.M      - Force same voxel-to-world to all images [read from file]
%
% optimise.pg.w     - Optimise subspace [true] or keep if fixed (false)
% optimise.z.z      - Optimise latent coordinates [true]
% optimise.z.A      - Optimise latent precision [true]
% optimise.q.q      - Optimise affine coordinates [true]
% optimise.q.A      - Optimise affine precision [true]
% optimise.r.r      - Optimise reidual fields [true]
% optimise.r.l      - Optimise residual precision [true]
% optimise.tpl.a    - Optimise template [true]
% optimise.mixreg.w - Optimise mixture weight [true]
% optimise.mixreg.a - Optimise mixture weight prior [true]
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
%                 iteration ['pg_result.mat']
% fnames.model  - Structure of filenames for all file arrays
% fnames.dat    - Structure of filenames for all file arrays
% ondisk.model  - Structure of logical for temporary array [default_ondisk]
% ondisk.dat    - "      "       "       "       "       "
% _________________________________________________________________________
%
% FORMAT [model, dat, opt] = pgra_model(opt, dat)
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
%  [Az0, nz0] --> (Az) --> (z) --> <v> <-- (r) <-- (lam) <-- [lam0, nl0]
%                                           ^
%                                           |
%                                          [L]
%
% Data part
% ---------
%
% <v> --> <phi> --> <psi> <-- <xi> <-- <R> <-- (q) <-- (Aq) <-- [Aq0, nq0]
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

    global_start = tic;
    fprintf([' ' repmat('-',1,78) ' \n']);
    str_started = sprintf('%20s || PGRA model started...', datestr(now));
    fprintf(['| ' str_started repmat(' ', 1, 80-3-length(str_started)) '|\n']);
    fprintf([' ' repmat('-',1,78) ' \n\n']);
    cleanupObj = onCleanup(@() goodbye(global_start));
        
    % -----------
    % Parse input
    % -----------
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
    setpath('pgra');
    
    % ---------------------------------------------------------------------
    %    Initialisation
    % ---------------------------------------------------------------------
    if ~cont
        
        % Default parameters + prepare structures & file arrays
        % -----------------------------------------------------------------
        [opt,dat,model] = pgra_model_input(input,opt);    % Read observed
        opt             = pgra_model_default(opt);        % Read options
        [opt,dat,model] = pgra_model_data(opt,dat,model); % Set arrays

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
        [~, opt.pg.ld] = spm_shoot_greens('kernel', double(opt.tpl.lat), double([opt.tpl.vs opt.pg.prm]), opt.pg.bnd);
        opt.pg.ld = opt.pg.ld(1);
        ker = spm_diffeo('kernel', double(opt.tpl.lat), double([opt.tpl.vs opt.pg.prm]));
        opt.pg.ker = [ker(1,1,1,1,1) ker(1,1,1,2,2) ker(1,1,1,3,3)];
        clear ker
        
        % Initialise all arrays (= model variables)
        % -----------------------------------------------------------------
        if opt.ui.verbose
            fprintf(['%10s | %10s | ' repmat('=',1,50) ' |\n'], 'EM', 'Init');
        end
        [dat, model] = pgra_model_init(dat, model, opt);
        
        % Some more stuff regarding processing
        % -----------------------------------------------------------------
        model.emit      = 1;
        model.q.active  = true;
        model.pg.active = true;
        model.r.active  = true;
        model.pg.ok     = 1;
        model.pg.ok2    = 0;
        model.pg.armijo = 1;
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
        pgra_plot_all(model, opt);
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
            elseif (opt.optimise.z.z || opt.optimise.pg.w) && ~model.pg.active
                model.pg.active = true;
                fprintf('%10s | %10s\n', 'Activate', 'PG');
            elseif opt.optimise.r.r && ~model.r.active
                model.r.active = true;
                fprintf('%10s | %10s\n', 'Activate', 'Residual');
            else
                fprintf('Converged :D\n');
                break
            end
        end
        model.emit = emit;
        
        % -----------------------------------------------------------------
        %    Affine
        % -----------------------------------------------------------------
        if model.q.active
        
            if opt.optimise.q.q
                % Update parameters
                % -----------------
                [dat, model] = pgra_batch('FitAffine', dat, model, opt);

                % -----------
                % Lower bound
                model = updateLowerBound(model);
                pgra_plot_all(model, opt);
                % -----------
            end
            
            if opt.q.Mr && opt.optimise.q.A
                % Update prior
                % ------------
                rind = opt.q.rind;
                model.q.A = spm_prob('Wishart', 'up', ...
                                     opt.N, 0, model.q.qq(rind,rind) + model.q.S(rind,rind), ...
                                     opt.q.A0, opt.q.n0);

                % -----------
                % Lower bound
                [dat, model]    = pgra_batch('LB', 'PrecisionQ', dat, model, opt);
                if opt.q.n0
                    model.lb.Aq.val = -spm_prob('Wishart', 'kl', ...
                                           model.Aq,         opt.nq0+opt.N, ...
                                           eye(numel(rind)), opt.nq0, ...
                                           'normal');
                end
                model = updateLowerBound(model);
                pgra_plot_all(model, opt);
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

                % Penalise previous failure
                % -------------------------
                if model.pg.ok < 0
                    model.pg.ok = model.pg.ok + 1;
                else

                    [dat, model] = pgra_batch('GradHessSubspace', dat, model, opt);

                    % Boundary conditions
                    % -------------------
                    spm_diffeo('boundary', opt.pg.bnd);

                    % Factor of the prior : ln p(z|W) + ln p(W)
                    % -------------------
                    reg = diag(model.mixreg.w(2) * (model.z.zz + model.z.S) + model.pg.n * eye(opt.pg.K));

                    model.pg.d = prepareOnDisk(model.pg.d, size(model.pg.w));
                    for k=1:opt.pg.K
                        % Search direction
                        % ----------------
                        model.pg.d(:,:,:,:,k) = -spm_diffeo('fmg', ...
                            single(model.pg.h(:,:,:,:,k)), single(model.pg.g(:,:,:,:,k)), ...
                            double([opt.tpl.vs reg(k) * opt.pg.prm 2 2]));
                        clear gw
                    end
                    model.pg.g = rmarray(model.pg.g);
                    model.pg.h = rmarray(model.pg.h);

                    [~, model, dat] = lsSubspace(model.pg.d, model, dat, opt);

                    if model.pg.ok > 0
                        model.pg.ok2 = 0;

                        % -----------
                        % Lower bound
                        model = updateLowerBound(model);
                        pgra_plot_all(model, opt);
                        % -----------
                    else
                        model.pg.ok2 = model.pg.ok2 - 1;
                        model.pg.ok  = model.pg.ok2;
                    end

                end % < penalise previous failure
                
            end

            % Latent coordinates
            % -------------------------------------------------------------
            if opt.optimise.z.z
                % Coordinates (variational update)
                % --------------------------------
                [dat, model] = pgra_batch('FitLatent', dat, model, opt);

                % -----------
                % Lower bound
                [dat, model] = pgra_batch('LB', 'Latent', dat, model, opt);
                model = updateLowerBound(model);
                pgra_plot_all(model, opt);
                % -----------
            end

            % Orthogonalisation
            % -------------------------------------------------------------
            if opt.optimise.pg.w && opt.optimise.z.z && opt.optimise.z.A
            
                % Orthogonalise
                % -------------
                if opt.ui.verbose, fprintf('%10s | %10s ', 'Ortho', ''); tic; end
                [U, iU] = orthogonalisationMatrix(spm_matcomp('LoadDiag', model.z.zz), spm_matcomp('LoadDiag', model.pg.ww));
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
                [dat, model] = pgra_batch('RotateSubspace', Q, iQ, dat, model, opt);
           
                % Latent precision
                % ----------------
                model.z.A = spm_prob('Wishart', 'up', ...
                                     opt.N, 0, model.z.zz + model.z.S, ...
                                     opt.z.A0, opt.z.n0);

                % -----------
                % Lower bound
                [dat, model] = pgra_batch('LB', 'Orthogonalise', dat, model, opt);
                model = updateLowerBound(model);
                pgra_plot_all(model, opt);
                % -----------
                

            % Latent precision
            % -------------------------------------------------------------
            elseif opt.optimise.z.A
                
                model.z.A = spm_prob('Wishart', 'up', ...
                                     opt.N, 0, model.z.zz + model.z.S, ...
                                     opt.z.A0, opt.z.n0);

                % -----------
                % Lower bound
                [dat, model] = pgra_batch('LB', 'PrecisionZ', dat, model, opt);
                model        = updateLowerBound(model);
                pgra_plot_all(model, opt);
                % -----------
            end
            
        end
        
        % -----------------------------------------------------------------
        %    Residual field
        % -----------------------------------------------------------------
        if model.r.active
            
            % Update subjects
            % ---------------
            if opt.optimise.r.r
                [dat, model] = pgra_batch('FitResidual', dat, model, opt);

                % -----------
                % Lower bound
                model = updateLowerBound(model);
                pgra_plot_all(model, opt);
                % -----------
            end
            
            % Update precision
            % ----------------
            if opt.optimise.r.l
                K = 3*prod(opt.tpl.lat);
                model.r.l = opt.r.n0/opt.r.l0 + (model.mixreg.w(1)/K)*(model.r.tr + model.r.reg);
                model.r.l = model.r.n/model.r.l;
                if opt.ui.verbose, fprintf('%10s | %10g\n', 'Lambda', model.r.l); end

                % -----------
                % Lower bound
                [dat, model]    = pgra_batch('LB', 'Lambda', dat, model, opt);
                model = updateLowerBound(model);
                pgra_plot_all(model, opt);
                % -----------
            end
            
        end
        
        % -----------------------------------------------------------------
        %    Template
        % -----------------------------------------------------------------
        if opt.optimise.tpl.a
            if opt.ui.verbose, fprintf('%10s | %10s ', 'Template', ''); tic; end
            if opt.tpl.cat
                model.tpl.a = updateMuML(opt.model, dat, ...
                    'lat', opt.tpl.lat, 'par', opt.split.par, ...
                    'debug', opt.ui.debug, 'output', model.tpl.a);
                model.tpl.gmu = templateGrad(model.tpl.a, ...
                    opt.tpl.itrp, opt.tpl.bnd,  ...
                    'debug',  opt.ui.debug, 'output', model.tpl.gmu);
                model.tpl.mu = reconstructProbaTemplate(model.tpl.a, ...
                    'par', opt.split.par, 'debug',  opt.ui.debug, ...
                    'output', model.tpl.mu);
            else
                model.tpl.mu = updateMuML(opt.model, dat, ...
                    'lat', opt.tpl.lat, 'par',    opt.split.par, ...
                    'debug',  opt.ui.debug, 'output', model.tpl.mu);
                model.tpl.gmu = templateGrad(model.tpl.mu, ...
                    opt.tpl.itrp, opt.tpl.bnd, ...
                    'debug',  opt.ui.debug, 'output', model.tpl.gmu);
            end
            if opt.ui.verbose, fprintf('| %6.3s\n', toc); end
            % -----------
            % Lower bound
            [dat, model] = pgra_batch('LB', 'Template', dat, model, opt);
            model = updateLowerBound(model);
            pgra_plot_all(model, opt);
            % -----------
        end
    end
    
end

% =========================================================================
function goodbye(global_start)
    
    global_end = toc(global_start);
    fprintf('\n');
    fprintf([' ' repmat('-',1,78) ' \n']);
    str_end_1 = sprintf('%s || PGRA model ended.', datestr(now));
    fprintf(['| ' str_end_1 repmat(' ', 1, 80-3-length(str_end_1)) '|\n']);
    str_end_2 = sprintf('%20s || ', 'Elapsed time');
    % Convert to units
    dur = duration(0,0,global_end);
    elapsed = floor(years(dur));
    dur = dur - years(elapsed(end));
    elapsed = [elapsed floor(days(dur))];
    dur = dur - days(elapsed(end));
    elapsed = [elapsed floor(hours(dur))];
    dur = dur - hours(elapsed(end));
    elapsed = [elapsed floor(minutes(dur))];
    dur = dur - minutes(elapsed(end));
    elapsed = [elapsed floor(seconds(dur))];
    units   = {'year' 'day' 'hour' 'minute' 'second'};
    for i=1:numel(elapsed)
        if elapsed(i) > 0
            str_end_2 = [str_end_2 sprintf('%d %s', elapsed(i), units{i})];
            if elapsed(i) > 1
                str_end_2 = [str_end_2 's'];
            end
            if sum(elapsed(i+1:end)) > 0
                str_end_2 = [str_end_2 ', '];
            end
        end
    end
    fprintf(['| ' str_end_2 repmat(' ', 1, 80-3-length(str_end_2)) '|\n']);
    fprintf([' ' repmat('-',1,78) ' \n\n']);
    diary off
    
end
