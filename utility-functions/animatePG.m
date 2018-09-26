function frames = animatePG(varargin)
% _________________________________________________________________________
%
% FORMAT animatePG(dat, model, opt, (fname), ...) > Animate subject
% FORMAT animatePG(model, opt, (fname), ...)      > Animate template
% FORMAT animatePG(model, (fname), ...)           > Manually built model
%
% -------------------------------------------------------------------------
%
% DATA/MODEL (obtained from pgra_model/pgva_model)
% ----------
% dat    - Data of one subject      (from shape_model)
% model  - Model                    (from shape_model or manually built)
% opt    - Option                   (from shape_model)
%
% Note: for manually built models, model should be a structure with fields:
%   A     - Latent precision matrix (or its diagonal)
%   w     - Principal subspace
%   prm   - Regularisation parameters for Geodesic shooting
%   a     - Template or Log-template
%   (vs)  - Voxel size [1 1 1]
%   (cat) - Is template categorical? [false]
%
% OPTIONAL ARGUMENTS
% ------------------
% fname  - Filename of the output video (.mp4, .avi, .gif) [do not write]
%
% KEYWORD ARGUMENTS
% -----------------
% pg       - Principal geodesic along which to shoot ([1], 2, ..., K)
% sigma    - Shooting range in terms of number of standard deviation [3]
% slice    - Slice direction (1 = sagittal, 2 = coronal, [3] = axial)
% x        - Slice location [sagittal: 1/3, axial/coronal: middle]
% n        - Number of frames generated [21]
% loop     - Repeat frames to form a loop? [false]
% duration - Total duration of the video in sec [5]
% middle   - Start from the middle frame (= no deformation) [false]
%
% OUTPUT
% ------
% frames - Structure array of Matlab frames. 
%          They can be given to the movie function.
% _________________________________________________________________________

    % ---------------------------------------------------------------------
    % Detect case
    if nargin < 3 || ~isstruct(varargin{3})
        shoot_subject = false;
        model         = varargin{1};
        if nargin < 2 || ~isstruct(varargin{2})
            opt = struct;
            varargin = varargin(2:end);
        else
            opt      = varargin{2};
            varargin = varargin(3:end);
        end
    else
        shoot_subject = true;
        dat           = varargin{1};
        model         = varargin{2};
        opt           = varargin{3};
        varargin      = varargin(4:end);
    end
    
    % ---------------------------------------------------------------------
    % Parse arguments
    p = inputParser;
    p.FunctionName = 'animatePG';
    p.addOptional('fname',     '',    @ischar);
    p.addParameter('pg',       1,     @(X) isscalar(X) && isnumeric(X));
    p.addParameter('sigma',    3,     @(X) isscalar(X) && isnumeric(X));
    p.addParameter('slice',    3,     @(X) isscalar(X) && isnumeric(X));
    p.addParameter('x',        nan,   @(X) isscalar(X) && isnumeric(X));
    p.addParameter('n',        21,    @(X) isscalar(X) && isnumeric(X));
    p.addParameter('loop',     false, @(X) isscalar(X) && (isnumeric(X)||islogical(X)));
    p.addParameter('duration', 5,     @(X) isscalar(X) && isnumeric(X));
    p.addParameter('middle',   false, @(X) isscalar(X) && (isnumeric(X)||islogical(X)));
    p.parse(varargin{:});
    pg       = p.Results.pg;
    fname    = p.Results.fname;
    sigma    = p.Results.sigma;
    slice    = p.Results.slice;
    x        = p.Results.x;
    n        = p.Results.n;
    loop     = p.Results.loop;
    duration = p.Results.duration;
    middle   = p.Results.middle;
    
    
    % ---------------------------------------------------------------------
    % Read model
    try
        % --- Model obtained from shape_model
        A   = model.z.A;
        w   = model.pg.w;
        vs  = opt.tpl.vs;
        prm = opt.pg.prm;
        iscat = opt.tpl.cat && checkarray(model.tpl.a);
        if iscat && checkarray(model.tpl.a)
            mu = model.tpl.a;
        else
            mu = model.tpl.mu;
        end
    catch
        try
            % --- "Manual" model
            A = model.A;
            if size(A,1) == 1 || size(A,2) == 1
                A = diag(A);
            end
            w = model.w;
            if isfield(model, 'vs')
                vs = model.vs;
            else
                vs = [1 1 1];
            end
            prm = model.prm;
            if isfield(model, 'cat')
                iscat = model.cat;
            else
                iscat = false;
            end
            mu = model.a;
        catch
            error('The model is not organised as expected');
        end
    end
    
    % ---------------------------------------------------------------------
    % Discretize z space
    cov  = inv(A);
    sd   = sqrt(cov(pg,pg));
    allz = linspace(-sigma*sd, sigma*sd, n);
    
    % ---------------------------------------------------------------------
    % Prepare output
    loops = numel(allz);
    frames(loops) = struct('cdata',[],'colormap',[]);
    i = 1;
    fprintf('Build and record frames\n');
    fig = figure;
    
    % ---------------------------------------------------------------------
    % If subject: start from its z coordinates
    if shoot_subject
        true_z = dat.z.z;
        mean   = reconstructVelocity('latent', true_z, 'subspace', w);
        r      = numeric(dat.v.v) - mean;
        clear mean
    end
    
    % ---------------------------------------------------------------------
    % Loop: build transformation and save frame
    for z=allz
        fprintf('z: %f\n', z);
        % --- Compute velocity
        if shoot_subject
            z1 = true_z;
            z1(pg) = z;
            v  = reconstructVelocity('latent', z1, 'subspace', w, 'residual', r);
        else
            v = w(:,:,:,:,pg) * z;
        end
        % --- Warp template
        iphi = exponentiateVelocity(v, 'iphi', 'vs', vs, 'prm', prm);
        mu1 = warp(iphi, mu);
        if iscat
            mu1 = reconstructProbaTemplate(mu1);
        end
        % --- Save frame
        h = imshow_categorical(mu1, vs, x, slice);
        frames(i) = getframe(h.Parent);
        i = i + 1;
    end
    close(fig);
    
    % ---------------------------------------------------------------------
    % Create final video
    fprintf('Save frames\n');
    
    if      loop &&  middle, iter = [ceil(n/2):n n:-1:1 1:ceil(n/2)];
    elseif  loop && ~middle, iter = [1:n n:-1:1];
    elseif ~loop &&  middle, iter = [ceil(n/2):n n:-1:ceil(n/2)];
    elseif ~loop && ~middle, iter = 1:n;
    end
    framerate = numel(iter)/duration;
    if endsWith(fname, 'mp4'), vidopt = {'MPEG-4'};
    else,                      vidopt = {};
    end
    if endsWith(fname, 'avi') || endsWith(fname, 'mp4')
        vid = VideoWriter(fname, vidopt{:});
        vid.FrameRate = framerate;
        open(vid);
    end
    first_frame = true;
    for i=iter
        fprintf('z: %f\n', allz(i));
        if endsWith(fname, 'gif')
            [imind, cm] = rgb2ind(frame2im(frames(i)), 256);
            if first_frame
                imwrite(imind, cm, fname, 'gif', 'Loopcount', inf, 'DelayTime', 1/framerate); 
                first_frame = false;
            else
                imwrite(imind, cm, fname, 'gif', 'WriteMode', 'append', 'DelayTime', 1/framerate); 
            end
        elseif endsWith(fname, 'avi') || endsWith(fname, 'mp4') 
            writeVideo(vid, frames(i));
        end
    end
    if endsWith(fname, 'avi') || endsWith(fname, 'mp4') 
        close(vid);
    end

end