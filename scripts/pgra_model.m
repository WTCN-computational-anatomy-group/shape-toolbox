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
%   q    = parameters of a rigid-body (or affine) transform, to align shapes
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
% pg.geod  - Additional geodesic prior on velocity fields [true]
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
% PROCESSING
% ----------
% optimise.pg.w  - Optimise subspace [true] or keep if fixed (false)
% optimise.z.z   - Optimise latent coordinates [true]
% optimise.z.A   - Optimise latent precision [true]
% optimise.q.q   - Optimise affine coordinates [true]
% optimise.q.A   - Optimise affine precision [true]
% optimise.r.r   - Optimise reidual fields [true]
% optimise.r.l   - Optimise residual precision [true]
% optimise.tpl.a - Optimise template [true]
% iter.em      - Maximum number of EM iterations [1000]
% iter.gn      - Maximum number of Gauss-Newton iterations [1]
% iter.ls      - Maximum number of line search iterations [6]
% iter.itg     - Number of integration steps for geodesic shooting [auto]
% lb.threshold - Convergence criterion (lower bound gain) [1e-5]
% split.loop   - How to split array processing: 'none'/'slice'/['subject']
% split.par    - Parallelise processing (number of workers): 0/n/[inf]
% split.batch  - Batch size for parallelisation [auto]
% ui.verbose   - Talk during processing [true]
% ui.debug     - Further debuging talk [false]
% ui.ftrack    - Figure object for the lower bound tracking [gcf]
% dist         - distributed processing. See `help distribute_default`.
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
    plotAll(model, opt);
    
    % ---------------------------------------------------------------------
    %    EM iterations
    % ---------------------------------------------------------------------
    for emit = model.emit:opt.iter.em
    
        % -----------------------------------------------------------------
        %    General tracking
        % -----------------------------------------------------------------
        
        if opt.ui.verbose
            fprintf(['%10s | %10d | ' repmat('=',1,50) ' |\n'], 'EM', emit);
        end
        
        if model.lb.lb.gain < opt.lb.threshold
            if opt.optimise.q.q && ~model.q.active
                model.q.active = true;
                fprintf('%10s | %10s\n', 'Activate', 'Affine');
            elseif (opt.optimise.pg.w || opt.optimise.z.z) && ~model.pg.active
                model.pg.active = true;
                fprintf('%10s | %10s\n', 'Activate', 'PG');
            elseif opt.optimise.r.r && ~model.r.active
                model.r.active = true;
                fprintf('%10s | %10s\n', 'Activate', 'Residual');
            else
                fprintf('%10s |\n', 'Converged :D');
                break
            end
        end
        model.emit = emit;
        
        % -----------------------------------------------------------------
        %    Affine
        % -----------------------------------------------------------------
        if opt.optimise.q.q && model.q.active
        
            % Update parameters
            % -----------------
            [dat, model] = pgra_batch('FitAffine', dat, model, opt);

            % -----------
            % Lower bound
            model.lb.q.list = [model.lb.q.list model.lb.q.val];
            model.lb.q.it   = [model.lb.q.it   emit];
            model = updateLowerBound(model);
            plotAll(model, opt);
            % -----------
            
            % Update prior
            % ------------
            if opt.optimise.q.A && opt.q.Mr
                rind = opt.q.rind;
                model.q.A = spm_prob('Wishart', 'up', ...
                                     opt.N, 0, model.q.qq(rind,rind) + model.q.S(rind,rind), ...
                                     opt.q.A0, opt.q.n0);

                % -----------
                % Lower bound
                [dat, model]    = pgra_batch('LB', 'PrecisionQ', dat, model, opt);
                model.lb.q.list = [model.lb.q.list model.lb.q.val];
                model.lb.q.it   = [model.lb.q.it   emit + 0.5];
                if opt.q.n0
                    model.lb.Aq.val = -spm_prob('Wishart', 'kl', ...
                                           model.Aq,         opt.nq0+opt.N, ...
                                           eye(numel(rind)), opt.nq0, ...
                                           'normal');
                    model.lb.Aq.list = [model.lb.Aq.list model.lb.Aq.val];
                    model.lb.Aq.it   = [model.lb.Aq.it   emit];
                end
                model.lb.m.list = [model.lb.m.list model.lb.m.val];
                model.lb.m.it   = [model.lb.m.it   emit];
                model = updateLowerBound(model);
                plotAll(model, opt);
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
                    reg = opt.pg.geod * (model.z.zz + model.z.S) + eye(opt.pg.K);
                    % /!\ Should I include Sz ?

                    model.pg.d = prepareOnDisk(model.pg.d, size(model.pg.w));
                    for k=1:opt.pg.K
                        % Gradient
                        % --------
                        lw = spm_diffeo('vel2mom', single(model.pg.w(:,:,:,:,k)), [opt.tpl.vs, opt.pg.prm]);
                        gw = single(model.pg.g(:,:,:,:,k)) + reg(k,k) * lw;
                        clear lw
                        
                        % Search direction
                        % ----------------
                        model.pg.d(:,:,:,:,k) = -spm_diffeo('fmg', ...
                            single(model.pg.h(:,:,:,:,k)), gw, ...
                            double([opt.tpl.vs reg(k,k) * opt.pg.prm 2 2]));
                        clear gw
                    end
                    model.pg.g = rmarray(model.pg.g);
                    model.pg.h = rmarray(model.pg.h);

                    [~, model, dat] = lsSubspace(model.pg.d, model, dat, opt);

                    if model.pg.ok > 0
                        model.pg.ok2 = 0;

                        % -----------
                        % Lower bound
                        model.lb.w.list = [model.lb.w.list model.lb.w.val];
                        model.lb.w.it   = [model.lb.w.it   emit];
                        model.lb.m.list = [model.lb.m.list model.lb.m.val];
                        model.lb.m.it   = [model.lb.m.it   emit + 1/5];
                        model.lb.z.list = [model.lb.z.list model.lb.z.val];
                        model.lb.z.it   = [model.lb.z.it   emit];
                        model           = updateLowerBound(model);
                        plotAll(model, opt);
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
                model.lb.m.list = [model.lb.m.list model.lb.m.val];
                model.lb.m.it   = [model.lb.m.it   emit + 2/5];
                model.lb.z.list = [model.lb.z.list model.lb.z.val];
                model.lb.z.it   = [model.lb.z.it   emit + 1/3];
                model           = updateLowerBound(model);
                plotAll(model, opt);
                % -----------
            end

            % Orthogonalisation
            % -------------------------------------------------------------
            if opt.optimise.z.z && opt.optimise.pg.w && opt.optimise.z.A
            
                % Orthogonalise
                % -------------
                if opt.ui.verbose, fprintf('%10s | %10s ', 'Ortho', ''); tic; end
                [U, iU] = orthogonalisationMatrix(spm_matcomp('LoadDiag', model.z.zz), spm_matcomp('LoadDiag', model.pg.ww));
                if opt.ui.verbose, fprintf('| %6.3fs\n', toc); end

                % Rescale
                % -------
                if opt.ui.verbose, fprintf('%10s | %10s ', 'Rescale', ''); tic; end
                [Q, iQ] = pgra_scale_pg(iU' * model.pg.ww * iU, ...
                                        U   * model.z.zz  * U', ...
                                        U   * model.z.S   * U', ...
                                        opt.z.A0, opt.z.n0, opt.N);
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
                [dat, model]    = pgra_batch('LB', 'Orthogonalise', dat, model, opt);
                model.lb.w.list = [model.lb.w.list model.lb.w.val];
                model.lb.w.it   = [model.lb.w.it   emit];
                if opt.z.n0
                    model.lb.Az.list = [model.lb.Az.list model.lb.Az.val];
                    model.lb.Az.it   = [model.lb.Az.it   emit];
                end
                model.lb.z.list = [model.lb.z.list model.lb.z.val];
                model.lb.z.it   = [model.lb.z.it   emit + 2/3];
                model = updateLowerBound(model);
                plotAll(model, opt);
                % -----------
                

            % Latent precision
            % -------------------------------------------------------------
            elseif opt.optimise.z.A
                
                model.z.A = spm_prob('Wishart', 'up', ...
                                     opt.N, 0, model.z.zz + model.z.S, ...
                                     opt.z.A0, opt.z.n0);

                % -----------
                % Lower bound
                [dat, model]    = pgra_batch('LB', 'PrecisionZ', dat, model, opt);
                if opt.z.n0
                    model.lb.Az.list = [model.lb.Az.list model.lb.Az.val];
                    model.lb.Az.it   = [model.lb.Az.it   emit];
                end
                model.lb.z.list = [model.lb.z.list model.lb.z.val];
                model.lb.z.it   = [model.lb.z.it   emit + 2/3];
                model = updateLowerBound(model);
                plotAll(model, opt);
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
                model.lb.m.list = [model.lb.m.list model.lb.m.val];
                model.lb.m.it   = [model.lb.m.it emit + 3/5];
                model.lb.r.list = [model.lb.r.list model.lb.r.val];
                model.lb.r.it   = [model.lb.r.it emit];
                model = updateLowerBound(model);
                plotAll(model, opt);
                % -----------
            end
            
            % Update precision
            % ----------------
            if opt.optimise.r.l
                K = 3*prod(opt.tpl.lat);
                model.r.l = opt.r.n0/opt.r.l0 + (1/K)*(model.r.tr + model.r.reg);
                model.r.l = (opt.N+opt.r.n0)/model.r.l;
                if opt.ui.verbose, fprintf('%10s | %10g\n', 'Lambda', model.r.l); end

                % -----------
                % Lower bound
                [dat, model]    = pgra_batch('LB', 'Lambda', dat, model, opt);
                model.lb.r.list = [model.lb.r.list model.lb.r.val];
                model.lb.r.it   = [model.lb.r.it emit + 1/4];
                if opt.r.n0
                    model.lb.l.list = [model.lb.l.list model.lb.l.val];
                    model.lb.l.it   = [model.lb.l.it   emit];
                end
                model = updateLowerBound(model);
                plotAll(model, opt);
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
                                         'lat',    opt.tpl.lat,   ...
                                         'par',    opt.split.par, ...
                                         'debug',  opt.ui.debug,  ...
                                         'output', model.tpl.a);
                model.tpl.gmu = templateGrad(model.tpl.a,  ...
                                             opt.tpl.itrp, ...
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
                                             opt.tpl.itrp, ...
                                             opt.tpl.bnd, ...
                                             'debug',  opt.ui.debug, ...
                                             'output', model.tpl.gmu);
            end
            if opt.ui.verbose, fprintf('| %6.3s\n', toc); end
            % -----------
            % Lower bound
            [dat, model] = pgra_batch('LB', 'Matching', dat, model, opt);
            model.lb.m.list = [model.lb.m.list model.lb.m.val];
            model.lb.m.it   = [model.lb.m.it   emit + 4/5];
            model = updateLowerBound(model);
            plotAll(model, opt);
            % -----------
        end
        
        % -----------------------------------------------------------------
        %    Lower bound gain
        % -----------------------------------------------------------------
        model = updateLowerBound(model, 'gain');
        
        % -----------------------------------------------------------------
        %    Save current state
        % -----------------------------------------------------------------
        if ~isempty(opt.fnames.result)
            createAllNifti(dat, model, opt);
            ftrack = opt.ui.ftrack;
            opt.ui.ftrack = nan;
            save(fullfile(opt.dir.model, opt.fnames.result), ...
                 'model', 'dat', 'opt');
            opt.ui.ftrack = ftrack;
        end
    end
end

% =========================================================================

function plotAll(model, opt)
% Plot PG + lower bound stuff
% This function is highly specific to this particular model. Not sure I can
% come up with a generic "plotModel" function, even though it could be
% nice.
    
    if opt.ui.verbose
        try
            figure(opt.ui.ftrack);
            clf(opt.ui.ftrack);
        catch
            figure(gcf);
            clf(gcf);
        end
        
        nw = 3;
        nh = 5;
        i  = 0;
        colors = ['b', 'g', 'r', 'c', 'm', 'k'];
        npg = nw;
        
        % --- Line 1
        
        % Template
        i = i + 1;
        subplot(nh, nw, i)
        tpl = catToColor(model.tpl.mu(:,:,ceil(size(model.tpl.mu,3)/2),:));
        dim = [size(tpl) 1 1];
        tpl = permute(reshape(tpl, [dim(1:2) dim(4)]), [2 1 3]);
        asp = 1./[opt.tpl.vs(2) opt.tpl.vs(1) 1];
        image(tpl(end:-1:1,:,:));
        daspect(asp);
        axis off
        npg = npg - 1;
        if size(model.tpl.mu, 3) > 1
            title('template (axial)')
            i = i + 1;
            subplot(nh, nw, i)
            tpl = catToColor(model.tpl.mu(:,ceil(size(model.tpl.mu,2)/2),:,:));
            dim = [size(tpl) 1 1];
            tpl = permute(reshape(tpl, [dim(1) dim(3) dim(4)]), [2 1 3]);
            asp = 1./[opt.tpl.vs(3) opt.tpl.vs(1) 1];
            image(tpl(end:-1:1,:,:));
            daspect(asp);
            axis off
            title('template (coronal)')
            npg = npg - 1;
        else
            title('template')
        end
        % PG
        for k=1:npg
            i = i + 1;
            subplot(nh,nw,i)
            pg = defToColor(model.pg.w(:,:,ceil(size(model.pg.w,3)/2),:,k));
            dim = [size(pg) 1 1];
            pg = permute(reshape(pg, [dim(1:2) dim(4)]), [2 1 3]);
            asp = 1./[opt.tpl.vs(2) opt.tpl.vs(1) 1];
            image(pg(end:-1:1,:,:));
            daspect(asp);
            axis off
            title(sprintf('PG %d', k))
        end
        
        % --- Line 2
        
        % Lower bound
        i = i + 1;
        subplot(nh,nw,i)
        plot([model.lb.lb.it model.lb.lb.curit], ...
             [model.lb.lb.list model.lb.lb.curlist], ...
             colors(mod(i, length(colors))+1))
        title('Lower bound')
        % Data likelihood
        i = i + 1;
        subplot(nh,nw,i)
        plot(model.lb.m.it, model.lb.m.list, ...
             colors(mod(i, length(colors))+1))
        title(model.lb.m.name)
        % PG prior
        i = i + 1;
        subplot(nh,nw,i)
        plot(model.lb.w.it, model.lb.w.list, ...
             colors(mod(i, length(colors))+1))
        title(model.lb.w.name)
        
        % --- Line 3
        
        if opt.q.Mr
            % KL affine
            i = i + 1;
            subplot(nh,nw,i)
            plot(model.lb.q.it, model.lb.q.list, ...
                 colors(mod(i, length(colors))+1))
            title(model.lb.q.name)
        else
            i = i + 1;
        end
        % KL residual
        i = i + 1;
        subplot(nh,nw,i)
        plot(model.lb.r.it, model.lb.r.list, ...
             colors(mod(i, length(colors))+1))
        title(model.lb.r.name)
        % KL latent
        i = i + 1;
        subplot(nh,nw,i)
        plot(model.lb.z.it, model.lb.z.list, ...
             colors(mod(i, length(colors))+1))
        title(model.lb.z.name)
        
        
        % --- Line 4
        
        if opt.q.Mr
            % KL affine precision
            i = i + 1;
            subplot(nh,nw,i)
            plot(model.lb.Aq.it, model.lb.Aq.list, ...
                 colors(mod(i, length(colors))+1))
            title(model.lb.Aq.name)
        else
            i = i + 1;
        end
        % KL residual precision
        i = i + 1;
        subplot(nh,nw,i)
        plot(model.lb.l.it, model.lb.l.list, ...
             colors(mod(i, length(colors))+1))
        title(model.lb.l.name)
        % KL latent precision
        i = i + 1;
        subplot(nh,nw,i)
        plot(model.lb.Az.it, model.lb.Az.list, ...
             colors(mod(i, length(colors))+1))
        title(model.lb.Az.name)
        
        % --- Line 5
        
        % WW
        i = i + 1;
        subplot(nh,nw,i)
        imagesc(model.pg.ww), colorbar;
        title('W''LW')
        % Sample covariance
        i = i + 1;
        subplot(nh,nw,i)
        imagesc(model.z.S + model.z.zz), colorbar;
        title('Sample covariance E[ZZ]')
        % Precision matrix
        i = i + 1;
        subplot(nh,nw,i)
        imagesc(model.z.A), colorbar;
        title('Latent precision E[A]')
        
        drawnow
    end

end