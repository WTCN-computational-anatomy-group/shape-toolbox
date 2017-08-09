classdef Velocity < handle & DiskWorker
% CLASS Velocity
% This class represents an initial velocity in terms of lattent coordinates  
% in the principal subspace plus a residual field.
% It allows to reconstruct a "true" velocity though V = WZ + s*R, where s
% is the ammount of "noise" in the population of velocities.
%
% MAP estimates of Z and R can be found by alternatingly maximizing
% P( X | W, Z, R ) P( Z | W, R ) P( Z ) and
% P( X | W, Z, R ) P( R | W, Z ) P( R ) 
% with respect to respectively Z and R with Gauss-Newton optimization.
% * P( Z | W, R ) = N( inv(W'LW) )
% * P( R | W, Z ) = N( inv(L) )
% * Matching term is one of normal (gaussian noise), laplace, binomial and 
% multinomial between a template and an image.
%
% # Public methods
% * Velocity('Key', value, ...)
% * update
%
% # Public properties
% ## Data
% > Z, Mu, Image, Vel, IPhi, ... (We can choose to make more stuff public)
% ## Options
% > MatchingTerm, RegParam, MaxIt, MaxGNIt.Z, MaxGNIt.R, Directory,
%   CanDistribute, MaxWorkers
%
% TODO:
% * I could also store the jacobians of direct and inverse transforms
% * Include affine registration
% * Take into account voxel2world mapping

    %% == Public properties ===============================================
    % --- Options
    properties
        MatchingTerm   = 'normal'       % [string]   Objective function class. One of normal, laplace, binomial, multinomial.
        RegParam       = [0.0001 0.001 0.2 0.05 0.2] % [5 double] Geodesic shooting parameters.
        SigmaR         = 1              % [double]   "Noise" captured by the residual field.
        MaxGNIt        = 1              % [int]      Max number of Gauss-Newton iterations.
        MaxLSIt        = 6              % [int]      Max number of line search iterations.
        Integration    = nan            % [int]      Number of steps for Euler-like integration
        Interpolation  = [1 1 1 1 1 1]  % [6 int]    Order/Boundary options for b-spline interpolation
        CanDistribute  = false          % [bool]     If true, Velocity's method may distribute loops for speed.
        MaxWorkers     = 4              % [int]      Maximum number of workers if parallelizing.
        Normal         = struct('s', 1) % [struct]   Additional options for 'normal' matching.
        Laplace        = struct('b', 1) % [struct]   Additional options for 'laplace' matching.
        Verbose        = false          % [bool]     Write additional information
        Graphic        = struct()       % [struct]   Set figures to see what's going on
        Debug          = false          % [bool]     Write super-additional information
    end
    % --- [Dependent] VoxelSize
    properties (Dependent, SetAccess = private)
        VoxelSize       % [3 double] Voxel size of the velocity lattice
    end
    methods
        function value = get.VoxelSize(obj)
            value = sqrt(sum(obj.Mu.mat(1:3,1:3).^2));
        end
    end
    % --- [Dependent] Directory
    properties (Dependent, GetAccess = private)
        Directory       % Change the directory of all nifti files
    end
    methods
        function set.Directory(obj, value)
            obj.disableListeners();
            array_names = fieldnames(obj.nii2dat);
            if isunix && ~startsWith(value, '/')
                value = [pwd value];
            elseif ispc && isempty(regexpi(value, '^[A-Z]+:[\\/]', 'start'))
                value = [pwd value];
            end
            for i=1:numel(array_names)
                [~, fname, ext] = fileparts(obj.(array_names{i}).dat.fname);
                obj.(array_names{i}).dat.fname = fullfile(value, [fname ext]);
            end
            obj.enableListeners();
        end
    end
    % --- Log-likelihood
    properties (GetObservable)
        ll              % [double] Complete model: ~= log p(F|Mu)
        llm             % [double] Matching term: log p(F|Mu, W, Z*, R*)
        llz             % [double] Z prior: log p(Z*|W)
        llr             % [double] R prior: log p(R*)
        lll             % [double] Laplace approximation
    end
    
    %% == Public niftis ===================================================
    % All "major" arrays that can be stored on disk are stored as nifti
    % images. We also add a few arrays which are not needed for the MAP
    % optimization but which make sense in terms of possible outputs. These
    % (transforms, warped images) will only be computed if they are asked
    % for. Public set/get access to these array consequently trigger
    % several listeners and is not super-efficient.
    % File arrays are directly access internally to avoid this call
    % overhead (see section "quick access to file_arrays").
    %--- (Possible) Input arrays ---
    properties (SetObservable, GetObservable)
        Z       % [nifti] Coordinates in the principal subspace
        W       % [nifti] Principal geodesic/component space
        R       % [nifti] Residual velocity field
        RegZ    % [nifti] Prior on Z (p(Z|W))
        Mu      % [nifti] Template
        GradMu  % [nifti] Spatial gradients of the template.
        LogMu   % [nifti] Log-probabilities (used with bi and multinomial matching terms)
        Image   % [nifti] Observed image
    end
    % --- Advanced setters ---
    % Setting a public array will trigger setNifti anyways, allowing to 
    % build/modify the appropriate nifti from an input which can be a
    % path/nifti/array...
    % These advance setters allow to pass advanced options to setNifti.
    % See > help setNifti < for available options.
    methods
        function setZ(obj, value, varargin)
        % FORMAT obj.setZ(value, 'OptKey', opt_value, ...)
        % See 'help setNifti' for a description of possible options.
            s = obj.disableListeners('Z');
            obj.Z = setNifti(value, obj.Z, varargin{:});
            obj.statusChanged('z');
            obj.enableListeners(s, 'Z'); 
        end
        function setW(obj, value, varargin)
        % FORMAT obj.setW(value, 'OptKey', opt_value, ...)
        % See 'help setNifti' for a description of possible options.
            s = obj.disableListeners('W');
            obj.W = setNifti(value, obj.W, varargin{:});
            obj.statusChanged('w');
            obj.enableListeners(s, 'W');
        end
        function setR(obj, value, varargin)
        % FORMAT obj.setR(value, 'OptKey', opt_value, ...)
        % See 'help setNifti' for a description of possible options.
            s = obj.disableListeners('R');
            obj.R = setNifti(value, obj.R, varargin{:});
            obj.statusChanged('r');
            obj.enableListeners(s, 'R');  
        end
        function setRegZ(obj, value, varargin)
        % FORMAT obj.setRegZ(value, 'OptKey', opt_value, ...)
        % See 'help setNifti' for a description of possible options.
            s = obj.disableListeners('RegZ');
            obj.RegZ = setNifti(value, obj.RegZ, varargin{:});
            obj.statusChanged('regz');
            obj.enableListeners(s, 'RegZ');
        end
        function setMu(obj, value, varargin)
        % FORMAT obj.setMu(value, 'OptKey', opt_value, ...)
        % See 'help setNifti' for a description of possible options.
            s = obj.disableListeners('Mu');
            obj.Mu = setNifti(value, obj.Mu, varargin{:});
            obj.statusChanged('mu');
            obj.enableListeners(s, 'Mu');
        end
        function setGradMu(obj, value, varargin)
        % FORMAT obj.setGradMu(value, 'OptKey', opt_value, ...)
        % See 'help setNifti' for a description of possible options.
            s = obj.disableListeners('GradMu');
            obj.GradMu = setNifti(value, obj.GradMu, varargin{:});
            obj.statusChanged('gmu');
            obj.enableListeners(s, 'GradMu');
        end
        function setImage(obj, value, varargin)
        % FORMAT obj.setImage(value, 'OptKey', opt_value, ...)
        % See 'help setNifti' for a description of possible options.
            s = obj.disableListeners('Image');
            obj.Image = setNifti(value, obj.Image, varargin{:});
            obj.statusChanged('f');
            obj.enableListeners(s, 'Image'); 
        end
        function setLogMu(obj, value, varargin)
        % FORMAT obj.setLogMu(value, 'OptKey', opt_value, ...)
        % See 'help setNifti' for a description of possible options.
            s = obj.disableListeners('setLogMu');
            obj.setLogMu = setNifti(value, obj.setLogMu, varargin{:});
            obj.statusChanged('a');
            obj.enableListeners(s, 'setLogMu');
        end
    end
    % --- Output arrays ---
    properties (GetObservable, SetAccess = private)
        Vel      % [nifti] Reconstructed velocity
        IPhi     % [nifti] Inverse transform (warps image to template ??)
        Phi      % [nifti] Direct transform (warps template to image ??)
        WMu      % [nifti] Warped template in image space
        WLogMu   % [nifti] Warped log-template in image space (used with bi and multinomial matching terms)
        WImage   % [nifti] Warped image in template space
        Res      % [nifti] Residuals between warped template and image
        PushedF  % [nifti] PushedImage in the template space
        PushedC  % [nifti] Number of image voxels pushed in each template voxel
        Hz       % [nifti] Hessian of LL w.r.t. Z (used for Laplace approximation)
        Hr       % [nifti] Hessian of LL w.r.t. R (used for Laplace approximation)
        Sz
        Sr
        Sv
        Swz
    end
    %% == Quick access to file_arrays =====================================
    properties (Dependent, Access = ?VelocityFriend)
        z        % [file_array] Coordinates in the principal subspace
        w        % [file_array] Principal geodesic/component space
        r        % [file_array] residual velocity field
        regz     % [file_array] Inverse covariance of the prior on Z (S^-1 in p(Z|W) = N(0, S))
        v        % [file_array] Reconstructed velocity
        iphi     % [file_array] Inverse transform (f = mu(iphi))
        phi      % [file_array] Direct transform (f(phi) = mu)
        mu       % [file_array] Template
        a        % [file_array] "Log-template" (encoding of the template in log/rotnull space)
        f        % [file_array] Observed image
        wmu      % [file_array] Warped template
        wa       % [file_array] Warped log-template
        wf       % [file_array] Warped observed image
        res      % [file_array] Residuals
        gmu      % [file_array] Spatial gradients of the template.
        pf       % [file_array] PushedImage in the template space
        pvox     % [file_array] Number of image voxels pushed in each template voxel
        hz       % [file_array] Hessian of LL w.r.t. Z (used for Laplace approximation)
        hr       % [file_array] Hessian of LL w.r.t. R (used for Laplace approximation)
        sz       % [file_array] Covariance matrix of Z's posterior (Laplace approximation)
        sr       % [file_array] Covariance matrix of R's posterior (Laplace approximation)
        swz      % [file_array] Covariance matrix of WZ's posterior (Laplace approximation)
        sv       % [file_array] Covariance matrix of Vel's posterior (Laplace approximation)
    end
    methods
        % --- file_array getters ---
        function val = get.gmu(obj),    val = obj.getDat('gmu');    end
        function val = get.pf(obj),     val = obj.getDat('pf');     end
        function val = get.pvox(obj),   val = obj.getDat('pvox');   end
        function val = get.hz(obj),     val = obj.getDat('hz');     end
        function val = get.hr(obj),     val = obj.getDat('hr');     end
        function val = get.sz(obj),     val = obj.getDat('sz');     end
        function val = get.sr(obj),     val = obj.getDat('sr');     end
        function val = get.sv(obj),     val = obj.getDat('sv');     end
        function val = get.swz(obj),    val = obj.getDat('swz');    end
        function val = get.z(obj),      val = obj.getDat('z');      end
        function val = get.w(obj),      val = obj.getDat('w');      end
        function val = get.r(obj),      val = obj.getDat('r');      end
        function val = get.regz(obj),   val = obj.getDat('regz');   end
        function val = get.v(obj),      val = obj.getDat('v');      end
        function val = get.iphi(obj),   val = obj.getDat('iphi');   end
        function val = get.phi(obj),    val = obj.getDat('phi');    end
        function val = get.mu(obj),     val = obj.getDat('mu');     end
        function val = get.f(obj),      val = obj.getDat('f');      end
        function val = get.wmu(obj),    val = obj.getDat('wmu');    end
        function val = get.wf(obj),     val = obj.getDat('wf');     end
        function val = get.res(obj),    val = obj.getDat('res');    end
        function val = get.a(obj)
            if any(strcmpi(obj.MatchingTerm, {'bernoulli', 'binomial', 'binary', 'categorical', 'multinomial'}))
                s = obj.disableListeners('LogMu');
                val = obj.LogMu.dat;
                obj.enableListeners(s, 'LogMu');
            else
                s = obj.disableListeners('Mu');
                val = obj.Mu.dat;
                obj.enableListeners(s, 'Mu');
            end
        end
        function val = get.wa(obj)
            if any(strcmpi(obj.MatchingTerm, {'bernoulli', 'binomial', 'binary', 'categorical', 'multinomial'}))
                s = obj.disableListeners('WLogMu');
                val = obj.WLogMu.dat;
                obj.enableListeners(s, 'WLogMu');
            else
                s = obj.disableListeners('WMu');
                val = obj.WMu.dat;
                obj.enableListeners(s, 'WMu');
            end
        end
        % --- file_array setters ---
        function set.gmu(obj, val),    obj.setDat('gmu', val);  end
        function set.pf(obj, val),     obj.setDat('pf', val);   end
        function set.pvox(obj, val),   obj.setDat('pvox', val); end
        function set.hz(obj, val),     obj.setDat('hz', val);   end
        function set.hr(obj, val),     obj.setDat('hr', val);   end
        function set.sz(obj, val),     obj.setDat('sz', val);   end
        function set.sr(obj, val),     obj.setDat('sr', val);   end
        function set.sv(obj, val),     obj.setDat('sv', val);   end
        function set.swz(obj, val),    obj.setDat('swz', val);  end
        function set.z(obj, val),      obj.setDat('z', val);    end
        function set.w(obj, val),      obj.setDat('w', val);    end
        function set.r(obj, val),      obj.setDat('r', val);    end
        function set.regz(obj, val),   obj.setDat('regz', val); end
        function set.v(obj, val),      obj.setDat('v', val);    end
        function set.iphi(obj, val),   obj.setDat('iphi', val); end
        function set.phi(obj, val),    obj.setDat('phi', val);  end
        function set.mu(obj, val),     obj.setDat('mu', val);   end
        function set.f(obj, val),      obj.setDat('f', val);    end
        function set.wmu(obj, val),    obj.setDat('wmu', val);  end
        function set.wf(obj, val),     obj.setDat('wf', val);   end
        function set.res(obj, val),    obj.setDat('res', val);  end
        function set.a(obj, val)
            if any(strcmpi(obj.MatchingTerm, {'bernoulli', 'binomial', 'binary', 'categorical', 'multinomial'}))
                s = obj.disableListeners('LogMu');
                obj.LogMu.dat = val;
                obj.enableListeners(s, 'LogMu');
            else
                s = obj.disableListeners('Mu');
                obj.Mu.dat = val;
                obj.enableListeners(s, 'Mu');
            end
        end
        function set.wa(obj, val)
            if any(strcmpi(obj.MatchingTerm, {'bernoulli', 'binomial', 'binary', 'categorical', 'multinomial'}))
                s = obj.disableListeners('WLogMu');
                obj.WLogMu.dat = val;
                obj.enableListeners(s, 'WLogMu');
            else
                s = obj.disableListeners('WMu');
                obj.WMu.dat = val;
                obj.enableListeners(s, 'WMu');
            end
        end
    end
    
    %% == Public methods ==================================================
    methods
        function obj = Velocity(varargin)
        % FORMAT obj = Velocity('Key', value, ...)
        %
        % ## Parameters
        % * 'MatchingTerm'  - Noise model Template>Image
        % * 'RegParam'      - Parameters for Geodesic shooting.
        % * 'SigmaR'        - "Noise" captured by the residual field.
        % * 'MaxGNIt'       - Maximum number of Gauss-Newton iterations.
        % * 'MaxLSIt'       - Maximum number of line search iterations.
        % * 'Interpolation' -
        % * 'Directory'     - Output directory for arrays stored on disk.
        % * 'CanDistribute' - Is the object authorized to distribute jobs?
        % * 'MaxWorkers'    - Maximum number of workers for parallelization.
        % * 'Normal.s'      - Noise variance for the gaussian model (sigma^2).
        % * 'Laplace.b'     - Noise variance for the laplace model.
        % * 'Verbose'       -
        % * 'Graphic'       -
        % * 'Debug'         -
        % * 'Directory'     - Default directory to store output files.
        %
        % ## Input arrays
        %
        % _The following arguments can take for value: an array, a Nifti 
        % or a path to a file._
        % * 'Z'     - Initial latent coordinates
        %             [default: 'random' for random initialization]
        % * 'W'     - Princial space (velocity/geodesic basis functions)
        %             [mandatory]
        % * 'R'     - Initial residual velocity field [default: null field]
        % * 'RegZ'  - Hessian of the prior term (p(Z|W)) w.r.t. Z.
        %            [mandatory for GN optimization]
        % * 'Mu'    - Template image [mandatory for GN optimization]
        % * 'Image' - Observed image [mandatory for GN optimization]
        %
        % _Additional arguments can be passed through (Z is used as
        % example). Default values are usually the same as in the input._
        % * 'Z.fname'   - File to use for working storage
        % * 'Z.dim'     - Array dimensions
        % * 'Z.dtype'   - Array data type
        % * 'Z.replace' - Completely replace the existing nifti object (true)
        %                 or just fill the old one (false).
        %                 [default: true]
            
            array_names = fieldnames(obj.nii2dat);
            
            % --- Define some default properties that can be post-applied
            popt = struct;
            for i=1:numel(array_names)
                popt.(array_names{i}).dtype = 'float32';
            end
            popt.RegZ.dtype     = 'float64';
            popt.PostCovZ.dtype = 'float64';

            % --- Convert ('ArrayName', value) to ('ArrayName.value', value)
            for i=1:numel(varargin)
                if ischar(varargin{i}) && any(strcmpi(varargin{i}, array_names))
                    varargin{i} = [varargin{i} '.value'];
                end
            end

            % --- Parse user arguments
            opt = struct('Directory', '');
            opt = parse_varargin(varargin, opt);
            
            % --- Post-define default properties
            if (isunix && ~startsWith(opt.Directory, '/')) || ...
               (ispc && isempty(regexpi(opt.Directory, '^[A-Z]+:[\\/]', 'start')))
                opt.Directory = fullfile(pwd, opt.Directory);
            end
            for i=1:numel(array_names)
                try
                    opt.(array_names{i}).value
                catch
                    try     opt.(array_names{i}).fname;
                    catch,  opt.(array_names{i}).fname = fullfile(opt.Directory, [array_names{i} '.nii']); end
                    try     opt.(array_names{i}).dtype;
                    catch,  opt.(array_names{i}).dtype = popt.(array_names{i}).dtype; end
                end
            end
            
            % --- Set properties
            function updateOption(obj, opt, option)
                if isfield(opt, option)
                    obj.(option) = opt.(option);
                end
            end
            updateOption(obj, opt, 'MatchingTerm');
            updateOption(obj, opt, 'MaxGNIt');
            updateOption(obj, opt, 'MaxLSIt');
            updateOption(obj, opt, 'Integration');
            updateOption(obj, opt, 'Interpolation');
            updateOption(obj, opt, 'CanDistribute');
            updateOption(obj, opt, 'MaxWorkers');
            updateOption(obj, opt, 'Normal');
            updateOption(obj, opt, 'Laplace');
            updateOption(obj, opt, 'SigmaR');
            updateOption(obj, opt, 'Verbose');
            updateOption(obj, opt, 'Graphic');
            updateOption(obj, opt, 'Debug');
            
            % --- Initialize arrays
            
            for i=1:numel(array_names) % < Public arrays
                obj.(array_names{i}) = setNifti([], opt.(array_names{i}));
                obj.utd.(obj.nii2dat.(array_names{i})) = isfield(opt.(array_names{i}), 'value');
            end
%             for i=1:numel(obj.private_nii) % < Private arrays
%                 basename = obj.nii2a.(obj.private_nii{i});
%                 fname = fullfile(opt.Directory, [basename '.nii']);
%                 obj.(obj.private_nii{i}) = setNifti([], 'fname', fname,  'dtype', 'float32' );
%                 obj.utd.(obj.nii2a.(obj.private_nii{i})) = false;
%             end
            
            % --- Try to guess some properties from file
            if ~isfield(opt, 'RegParam')
                trydesc = {'Z', 'W', 'R', 'Vel'};
                pattern = '(?<type>[\w\s]*) \((?<param>[^\(\)]*)\)';
                for i=1:numel(trydesc)
                    s = regexp(obj.(trydesc{i}).descrip, pattern, 'names');
                    if ~isempty(s)
                        par = str2num(s.param);
                        if length(par) == 5
                            obj.RegParam = par;
                            break
                        end
                    end
                end
            end
            
            % --- Listeners
            obj.setAllListeners()
        end
        
        % -- Externaly defined
        varargout = update(obj, varargin)
        varargout = updateZ(obj, varargin)
        varargout = updateR(obj, varargin)
        varargout = initZ(obj, varargin)
        varargout = initR(obj, varargin)
    end
    methods (Access = protected)
        % -- Externaly defined
        varargout = setDat(obj, varargin)
        varargout = getDat(obj, varargin)
        varargout = checkarray(obj, varargin)
        varargout = reconstructVelocity(obj, varargin)
        varargout = exponentiateVelocity(obj, varargin)
        varargout = pushImage(obj, varargin)
        varargout = warpTemplate(obj, varargin)
        varargout = warpImage(obj, varargin)
        varargout = reconstructTemplate(obj, varargin)
        varargout = reconstructWarpedTemplate(obj, varargin)
        varargout = computeTemplateGrad(obj, varargin)
        varargout = computeRegZ(obj, varargin)
        varargout = computeGradHessZ(obj, varargin)
        varargout = gradHessMatchingVel(obj, varargin)
        varargout = gradHessMatchingZ(obj, varargin)
        varargout = gradHessMatchingR(obj, varargin)
        varargout = gradHessPriorZ(obj, varargin)
        % varargout = gradHessRegZ(obj, varargin) % Additional reg
        varargout = logLikelihood(obj, varargin)         % TODO ?
        varargout = logLikelihoodMatching(obj, varargin) % log p(F | Z, R, W)
        varargout = logLikelihoodPriorZ(obj, varargin)   % log p(Z | W)
        varargout = logLikelihoodPriorR(obj, varargin)   % log p(R)
        varargout = gaussNewtonZ(obj, varargin)
        varargout = gaussNewtonR(obj, varargin)
        varargout = lineSearchZ(obj, varargin)
        varargout = lineSearchR(obj, varargin)
    end
    %% == Define DiskWorker properties ====================================
    properties (Constant, Access = protected)
        dependencies = struct(...
            'z',    {{}}, ...
            'w',    {{}}, ...
            'r',    {{}}, ...
            'regz', {{'w', 'RegParam'}}, ...
            'mu',   {{'a'}}, ...
            'f',    {{}}, ...
            'v',    {{'z', 'w', 'r'}}, ...
            'iphi', {{'v', 'RegParam', 'Integration'}}, ...
            'phi',  {{'v', 'RegParam', 'Integration'}}, ...
            'wmu',  {{'wa'}}, ...
            'wa',   {{'a', 'iphi', 'Interpolation'}}, ...
            'res',  {{'f', 'wmu', 'MatchingTerm', 'Normal', 'Laplace'}}, ...
            'wf',   {{'f', 'phi', 'Interpolation'}}, ...
            'gmu',  {{'a', 'Interpolation'}}, ...
            'pf',   {{'f', 'iphi'}}, ...
            'pvox', {{'f', 'iphi'}}, ...
            'hz',   {{'w', 'pf', 'pvox', 'gmu', 'mu', 'a', 'MatchingTerm', 'Normal', 'Laplace'}}, ...
            'hr',   {{'w', 'pf', 'pvox', 'gmu', 'mu', 'a', 'MatchingTerm', 'Normal', 'Laplace'}}, ...
            'sz',   {{'hz'}}, ...
            'sr',   {{'hr', 'RegParam'}}, ...
            'swz',  {{'sz', 'w'}}, ...
            'sv',   {{'sr', 'swz'}}, ...
            'llm',  {{'res', 'f', 'wmu', 'pf', 'pvox', 'MatchingTerm', 'Normal', 'Laplace'}}, ...
            'llz',  {{'z', 'regz'}}, ...
            'llr',  {{'r', 'RegParam'}}, ...
            'lll',  {{'hz', 'hr'}}, ...
            'll',   {{'llm', 'llz', 'llr', 'lll'}} ...
        )
        
        nii2dat = struct(...
            'Z',        'z',    ...
            'W',        'w',    ...
            'R',        'r',    ...
            'RegZ',     'regz', ...
            'Mu',       'mu',   ...
            'LogMu',    'a',    ...
            'Image',    'f',    ...
            'Vel',      'v',    ...
            'IPhi',     'iphi', ...
            'Phi',      'phi',  ...
            'WMu',      'wmu',  ...
            'WLogMu',   'wa',   ...
            'Res',      'res',  ...
            'WImage',   'wf',   ...
            'GradMu',   'gmu',  ...
            'PushedF',  'pf',   ...
            'PushedC',  'pvox', ...
            'Hz',       'hz',   ...
            'Hr',       'hr',   ...
            'Sz',       'sz',   ...
            'Sr',       'sr',   ...
            'Swz',      'swz',  ...
            'Sv',       'sv'    )
        
        updaters = struct(...
            'z',    @(X){}, ...
            'w',    @(X){}, ...
            'r',    @(X){}, ...
            'regz', @computeRegZ, ...
            'mu',   @reconstructTemplate, ...
            'a',    @(X){}, ...
            'f',    @(X){}, ...
            'v',    @reconstructVelocity, ...
            'iphi', @(X)exponentiateVelocity(X,'iphi'), ...
            'phi',  @(X)exponentiateVelocity(X,'phi'), ...
            'wmu',  @reconstructWarpedTemplate, ...
            'wa',   @warpTemplate, ...
            'res',  @residuals, ...
            'wf',   @warpImage, ...
            'gmu',  @computeTemplateGrad, ...
            'pf',   @pushImage, ...
            'pvox', @pushImage, ...
            'hz',   @(X){}, ...
            'hr',   @(X){}, ...
            'sz',   @uncertaintyZ, ...
            'sr',   @uncertaintyR, ...
            'swz',  @uncertaintyWZ, ...
            'sv',   @uncertaintyVel, ...
            'llm',  @logLikelihoodMatching, ...
            'llz',  @logLikelihoodPriorZ, ...
            'llr',  @logLikelihoodPriorR, ...
            'lll',  @logLikelihoodLaplace, ...
            'll',   @logLikelihood ...
        )
    
        dat2nii = DiskWorker.invertTranslation(Velocity.nii2dat)
        downdep = DiskWorker.invertDependencies(Velocity.dependencies)
    end
end