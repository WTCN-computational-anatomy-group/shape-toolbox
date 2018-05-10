function [model, dat, opt] = shape_model(varargin)
% _________________________________________________________________________
%
%                         Shape model (2018)
% _________________________________________________________________________
%
% FORMAT [model, dat, opt] = shape_model(input, (opt))
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

    cleanupObj = hello('Shape model');
    
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
    setpath('shape_model');
    
    
    % =====================================================================
    %    Initialisation
    % =====================================================================
    if ~cont
        
        % Default parameters + prepare structures & file arrays
        % -----------------------------------------------------------------
        [opt,dat,model] = shape_input(input,opt);    % Read observed
        opt             = shape_default(opt);        % Read options
        [opt,dat,model] = shape_data(opt,dat,model); % Set arrays

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
        opt.pg.LogDetL  = ldapprox(opt.pg.prm,  'vs', opt.tpl.vs, 'dim', [opt.tpl.lat 3], 'type', 'diffeo');
        opt.tpl.LogDetL = ldapprox(opt.tpl.prm, 'vs', opt.tpl.vs, 'dim', opt.tpl.lat,     'type', 'field');
        
        % Initialise all arrays (= model variables)
        % -----------------------------------------------------------------
        [dat, model] = shape_init(dat, model, opt);
        
        % Some more stuff regarding processing
        % -----------------------------------------------------------------
        model.emit = 0;
    end
    
    
    % =====================================================================
    %    Processing
    % =====================================================================
    switch lower(opt.par.subjects.mode)
        case 'qsub.tree'
            
            % -------------------------------------------------------------
            %    MODE :: CLUSTER TREE
            % -------------------------------------------------------------
            % [TODO]
            
        case {'qsub.flow', 'parfor', 'for'}
            
            % -------------------------------------------------------------
            %    MODE :: CLUSTER FLOW / PARFOR / FOR
            % -------------------------------------------------------------
            % We use the distribute toolbox to choose between:
            % qsub.flow = multiple matlab instances (cluster or workstation)
            % parfor    = Matlab's parfor 
            % for       = no parallelisation
            
            for emit=model.emit:opt.iter.em
                
                model = updateLowerBound(model, 'gain');
                
                if opt.lb.exact
                    
                    [dat, model] = shape_process(dat, model, opt);
                    
                else % Most efficient 
                
                    % Subject-specific processing
                    % -----------------------------------------------------
                    [opt.par.subjects, dat] = distribute(opt.par.subjects, ...
                        'shape_process_subject', 'iter', dat, model, opt);

                    % Population-specific processing
                    % -----------------------------------------------------
                    model = shape_process_pop(dat, model, opt);
                    
                end
                
                % Check convergence
                % ---------------------------------------------------------
                if model.converged
                    break
                end
                
            end
    end
       
end