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
%   q    = parameters of a rigid-body (or affine) transform, to align shapes
%   a/mu = Template (i.e., mean shape)
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
% pg.K     - Number of principal geodesics [32]
% pg.prm   - Parameters of the geodesic operator [1e-4 1e-3 0.2 0.05 0.2]
% pg.bnd   - Boundary conditions for the geodesic operator [1 = circulant]
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
% PROCESSING
% ----------
% optimise.pg.w  - Optimise subspace [true] or keep if fixed (false)
% optimise.z.A   - Optimise latent precision [true]
% optimise.q.A   - Optimise affine precision [true]
% optimise.v.l   - Optimise residual precision [true]
% optimise.tpl.a - Optimise template [true]
% iter.em      - Maximum number of EM iterations [1000]
% iter.gn      - Maximum number of Gauss-Newton iterations [1]
% iter.ls      - Maximum number of line search iterations [6]
% iter.itg     - Number of integration steps for geodesic shooting [auto]
% lb.threshold - Convergence criterion (lower bound gain) [1e-5]
% split.loop   - How to split array processing: 'none'/'slice'/['subject']
% split.par    - Parallelise processing (number of workers): 0/n/[inf]
% split.batch  - Batch size for parallelisation [auto]
% dist         - [TODO] distribute processing on a cluster
% ui.verbose   - Talk during processing [true]
% ui.debug     - Further debuging talk [false]
% ui.ftrack    - Figure object for the lower bound tracking [gcf]
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
    path = fileparts(which('pgva_model'));
    addpath(fullfile(path, 'pgva'));
    
    if ~cont
        % -----------------------------------------------------------------
        %    Default parameters + prepare structures & file arrays
        % -----------------------------------------------------------------
        [opt,dat,model] = pgva_model_input(input,opt);    % Read observed
        opt             = pgva_model_default(opt);        % Read options
        [opt,dat,model] = pgva_model_data(opt,dat,model); % Set arrays

        % -----------------------------------------------------------------
        %    Post-set parameters
        % -----------------------------------------------------------------
        spm_diffeo('boundary', opt.pg.bnd);
        [~, opt.pg.ld] = spm_shoot_greens('kernel', double(opt.tpl.lat), double([opt.tpl.vs opt.pg.prm]), opt.pg.bnd);
        opt.pg.ld = opt.pg.ld(1);
        
        % -----------------------------------------------------------------
        %    Initialise all variables
        % -----------------------------------------------------------------
        if opt.ui.verbose
            fprintf(['%10s | %10s | ' repmat('=',1,50) ' |\n'], 'EM', 'Init');
        end
        [dat, model] = pgva_model_init(dat, model, opt);
    end
    plotAll(model, opt);
    
    model.q.active  = false;
    model.v.active  = true;
    model.pg.active = true;
    model = updateLowerBound(model, 'gain');
    
    % ---------------------------------------------------------------------
    %    EM iterations
    % ---------------------------------------------------------------------
    for emit = 1:opt.iter.em
    
        % -----------------------------------------------------------------
        %    General tracking
        % -----------------------------------------------------------------
        
        if opt.ui.verbose
            fprintf(['%10s | %10d | ' repmat('=',1,50) ' |\n'], 'EM', emit);
        end
        
        if model.lb.lb.gain < opt.lb.threshold
            if ~model.q.active
                model.q.active = true;
                fprintf('%10s | %10s\n', 'Activate', 'Affine');
            elseif ~model.v.active
                model.v.active = true;
                fprintf('%10s | %10s\n', 'Activate', 'Velocity');
            elseif ~model.pg.active
                model.pg.active = true;
                fprintf('%10s | %10s\n', 'Activate', 'PG');
            else
                fprintf('%10s |\n', 'Converged :D');
                break
            end
        end
        model.emit = emit;
        
        % -----------------------------------------------------------------
        %    Affine
        % -----------------------------------------------------------------
        
        if model.q.active && opt.f.N
        
            % Update parameters
            % -----------------
            [dat, model] = pgva_batch('FitAffine', dat, model, opt);

            % Update prior
            % ------------
            rind = opt.q.rind;
            model.q.A = spm_prob('Wishart', 'up', ...
                                 opt.f.N, 0, model.q.qq(rind,rind) + model.q.S(rind,rind), ...
                                 opt.q.A0, opt.q.n0);

            % -----------
            % Lower bound
            model.lb.q.list = [model.lb.q.list model.lb.q.val];
            model.lb.q.it   = [model.lb.q.it   emit];
            if opt.q.n0
                model.lb.Aq.val = -spm_prob('Wishart', 'kl', ...
                                       model.Aq,         opt.nq0+opt.f.N, ...
                                       eye(numel(rind)), opt.nq0, ...
                                       'normal');
                model.lb.Aq.list = [model.lb.Aq.list model.lb.Aq.val];
                model.lb.Aq.it   = [model.lb.Aq.it   emit];
            end
            model = updateLowerBound(model);
            plotAll(model, opt);
            % -----------
            
        end
        
        % -----------------------------------------------------------------
        %    Velocity field
        % -----------------------------------------------------------------
        
        if model.v.active
        
            % Update subjects
            % ---------------
            [dat, model] = pgva_batch('FitVelocity', dat, model, opt);

            % -----------
            % Lower bound
            if opt.f.N
                model.lb.m.list = [model.lb.m.list model.lb.m.val];
                model.lb.m.it   = [model.lb.m.it emit];
                model.lb.v1.list = [model.lb.v1.list model.lb.v1.val];
                model.lb.v1.it   = [model.lb.v1.it emit];
            end
            if opt.v.N
                model.lb.v2.list = [model.lb.v2.list model.lb.v2.val];
                model.lb.v2.it   = [model.lb.v2.it emit];
            end
            model = updateLowerBound(model);
            plotAll(model, opt);
            % -----------
            
            % Update precision
            % ----------------
            K = 3*prod(opt.tpl.lat);
            model.v.l = opt.v.n0/opt.v.l0 + (1/K)*(model.v.uncty + model.v.tr + model.v.reg);
            model.v.l = (opt.f.N+opt.v.N+opt.v.n0)/model.v.l;
            if opt.ui.verbose, fprintf('%10s | %10g\n', 'Lambda', model.v.l); end
            
            % -----------
            % Lower bound
            [dat, model]    = pgva_batch('LB', 'Lambda', dat, model, opt);
            if opt.f.N
                model.lb.v1.list = [model.lb.v1.list model.lb.v1.val];
                model.lb.v1.it   = [model.lb.v1.it emit+1/4];
            end
            if opt.v.N
                model.lb.v2.list = [model.lb.v2.list model.lb.v2.val];
                model.lb.v2.it   = [model.lb.v2.it emit+1/4];
            end
            if opt.v.n0
                model.lb.l.list = [model.lb.l.list model.lb.l.val];
                model.lb.l.it   = [model.lb.l.it   emit];
            end
            model = updateLowerBound(model);
            plotAll(model, opt);
            % -----------
            
        end
        
        if model.pg.active
        
            % -------------------------------------------------------------
            %    Principal subspace
            % -------------------------------------------------------------

            if opt.ui.verbose, fprintf('%10s | %10s ', 'PG', ''); tic; end
            M = model.z.S + model.z.zz + eye(opt.pg.K) / model.v.l;
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
            [dat, model]    = pgva_batch('LB', 'Subspace', dat, model, opt);
            model.lb.w.list = [model.lb.w.list model.lb.w.val];
            model.lb.w.it   = [model.lb.w.it   emit];
            if opt.f.N
                model.lb.v1.list = [model.lb.v1.list model.lb.v1.val];
                model.lb.v1.it   = [model.lb.v1.it emit+2/4];
            end
            if opt.v.N
                model.lb.v2.list = [model.lb.v2.list model.lb.v2.val];
                model.lb.v2.it   = [model.lb.v2.it emit+2/4];
            end
            model = updateLowerBound(model);
            plotAll(model, opt);
            % -----------

            % -------------------------------------------------------------
            %    Latent coordinates
            % -------------------------------------------------------------

            % Coordinates (variational update)
            % --------------------------------
            [dat, model] = pgva_batch('FitLatent', dat, model, opt);

            % -----------
            % Lower bound
            [dat, model]    = pgva_batch('LB', 'Latent', dat, model, opt);
            model.lb.z.list = [model.lb.z.list model.lb.z.val];
            model.lb.z.it   = [model.lb.z.it   emit];
            if opt.f.N
                model.lb.v1.list = [model.lb.v1.list model.lb.v1.val];
                model.lb.v1.it   = [model.lb.v1.it emit+3/4];
            end
            if opt.v.N
                model.lb.v2.list = [model.lb.v2.list model.lb.v2.val];
                model.lb.v2.it   = [model.lb.v2.it emit+3/4];
            end
            model = updateLowerBound(model);
            plotAll(model, opt);
            % -----------

            % Orthogonalise
            % -------------
            if opt.ui.verbose, fprintf('%10s | %10s ', 'Ortho', ''); tic; end
            [U, iU] = orthogonalisationMatrix(spm_matcomp('LoadDiag', model.z.zz), spm_matcomp('LoadDiag', model.pg.ww));
            if opt.ui.verbose, fprintf('| %6.3fs\n', toc); end

            % Rescale
            % -------
            if opt.ui.verbose, fprintf('%10s | %10s ', 'Rescale', ''); tic; end
            [Q, iQ] = pgva_scale_pg(iU' * model.pg.ww * iU, ...
                                U   * model.z.zz  * U', ...
                                U   * model.z.S   * U', ...
                                model.v.l, opt.z.A0, opt.z.n0, opt.f.N+opt.v.N);
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
%                              
%             % -----------
%             % Lower bound
%             [dat, model] = pgva_batch('LB', 'PrecisionZ', dat, model, opt);
%             if opt.z.n0
%                 model.lb.Az.list = [model.lb.Az.list model.lb.Az.val];
%                 model.lb.Az.it   = [model.lb.Az.it   emit];
%             end
%             model.lb.z.list = [model.lb.z.list model.lb.z.val];
%             model.lb.z.it   = [model.lb.z.it   emit];
%             model = updateLowerBound(model);
%             plotAll(model, opt);
%             % -----------
                             
            % -----------
            % Lower bound
            [dat, model]    = pgva_batch('LB', 'PrecisionZ', dat, model, opt);
            [dat, model]    = pgva_batch('LB', 'Orthogonalise', dat, model, opt);
            model.lb.w.list = [model.lb.w.list model.lb.w.val];
            model.lb.w.it   = [model.lb.w.it   emit];
            if opt.z.n0
                model.lb.Az.list = [model.lb.Az.list model.lb.Az.val];
                model.lb.Az.it   = [model.lb.Az.it   emit];
            end
            if opt.f.N
                model.lb.v1.list = [model.lb.v1.list model.lb.v1.val];
                model.lb.v1.it   = [model.lb.v1.it emit+3/4];
            end
            if opt.v.N
                model.lb.v2.list = [model.lb.v2.list model.lb.v2.val];
                model.lb.v2.it   = [model.lb.v2.it emit+3/4];
            end
            model.lb.z.list = [model.lb.z.list model.lb.z.val];
            model.lb.z.it   = [model.lb.z.it   emit];
            model = updateLowerBound(model);
            plotAll(model, opt);
            % -----------
            
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
            [dat, model] = pgva_batch('LB', 'Matching', dat, model, opt);
            model.lb.m.list = [model.lb.m.list model.lb.m.val];
            model.lb.m.it   = [model.lb.m.it   emit];
            model = updateLowerBound(model);
            plotAll(model, opt);
            % -----------
        end
        
        % -----------------------------------------------------------------
        %    Save current state
        % -----------------------------------------------------------------
        model = updateLowerBound(model, 'gain');
        if ~isempty(opt.fnames.result)
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
        figure(opt.ui.ftrack);
        clf(opt.ui.ftrack);
        
        nw = 3;
        nh = 5;
        i  = 0;
        colors = ['b', 'g', 'r', 'c', 'm', 'k'];
        
        % --- Line 1
        
        % Template
        if opt.f.N
            i = i + 1;
            subplot(nh, nw, i)
            tpl = catToColor(model.tpl.mu(:,:,ceil(size(model.tpl.mu,3)/2),:));
            dim = [size(tpl) 1 1];
            image(reshape(tpl, [dim(1:2) dim(4)]));
            daspect(1./opt.tpl.vs);
            axis off
            title('template')
            npg = nw - 1;
        else
            npg = nw;
        end 
        % PG
        for k=1:npg
            i = i + 1;
            subplot(nh,nw,i)
            pg = defToColor(model.pg.w(:,:,ceil(size(model.pg.w,3)/2),:,k));
            dim = [size(pg) 1 1];
            image(reshape(pg, [dim(1:2) dim(4)]));
            daspect(1./opt.tpl.vs);
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
        if opt.f.N
            % Data likelihood
            i = i + 1;
            subplot(nh,nw,i)
            plot(model.lb.m.it, model.lb.m.list, ...
                 colors(mod(i, length(colors))+1))
            title(model.lb.m.name)
        else
            i = i +1;
        end
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
        elseif opt.v.N
            % LL residual
            i = i + 1;
            subplot(nh,nw,i)
            plot(model.lb.v2.it, model.lb.v2.list, ...
                 colors(mod(i, length(colors))+1))
            title(model.lb.v2.name)
        else
            i = i + 1;
        end
        if opt.f.N
            % KL residual
            i = i + 1;
            subplot(nh,nw,i)
            plot(model.lb.v1.it, model.lb.v1.list, ...
                 colors(mod(i, length(colors))+1))
            title(model.lb.v1.name)
        else
            i = i + 1;
        end
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

% =========================================================================
