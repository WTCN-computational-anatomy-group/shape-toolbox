function frames = animatePG(varargin)
% _________________________________________________________________________
%
% FORMAT (frames) = animatePG((dat), model, opt, (pg), (fname), ...)
%
% -------------------------------------------------------------------------
%
% DATA/MODEL (obtained from pgra_model/pgva_model)
% ----------
% dat    - Data of one subject
% model  - Model
% opt    - Option
%
% OPTIONAL
% --------
% pg     - Principal geodesic along which to shoot ([1], 2, ..., K)
% fname  - Filename of the output video (.mp4, .avi, .gif) [do not write]
%
% KEYWORD ARGUMENTS
% -----------------
% sigma - Shooting range in terms of number of standard deviation [3]
% slice  - Slice direction (1 = sagittal, 2 = coronal, [3] = axial)
% x      - Slice location [sagittal: 1/3, axial/coronal: middle]
% n      - Number of frames generated [21]
%
% OUTPUT
% ------
% frames - Structure array of Matlab frames. 
%          They can be given to the movie function.
% _________________________________________________________________________

    % Detect case
    % -----------
    if nargin < 3 || ~isstruct(varargin{3})
        shoot_subject = false;
        model = varargin{1};
        opt   = varargin{2};
        varargin = varargin(3:end);
    else
        shoot_subject = true;
        dat   = varargin{1};
        model = varargin{2};
        opt   = varargin{3};
        varargin = varargin(4:end);
    end

    % Parse arguments
    % ---------------
    p = inputParser;
    p.FunctionName = 'animatePG';
    p.addOptional('pg',    1,    @(X) isscalar(X) && isnumeric(X));
    p.addOptional('fname', '',   @ischar);
    p.addParameter('sigma', 3,   @(X) isscalar(X) && isnumeric(X));
    p.addParameter('slice', 3,   @(X) isscalar(X) && isnumeric(X));
    p.addParameter('x',     nan, @(X) isscalar(X) && isnumeric(X));
    p.addParameter('n',     21,  @(X) isscalar(X) && isnumeric(X));
    p.parse(varargin{:});
    pg    = p.Results.pg;
    fname = p.Results.fname;
    sigma = p.Results.sigma;
    slice = p.Results.slice;
    x     = p.Results.x;
    n     = p.Results.n;
    
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
    
    cov  = inv(A);
    sd   = sqrt(cov(pg,pg));
    allz = linspace(-sigma*sd, sigma*sd, n);
    
    loops = numel(allz);
    frames(loops) = struct('cdata',[],'colormap',[]);
    i = 1;
    fprintf('Build and record frames\n');
    fig = figure;
    if shoot_subject
        true_z   = dat.z.z;
        mean = reconstructVelocity('latent', true_z, 'subspace', w);
        r    = numeric(dat.v.v) - mean;
        clear mean
    end
    for z=allz
        fprintf('z: %f\n', z);
        if shoot_subject
            z1 = true_z;
            z1(pg) = z;
            v  = reconstructVelocity('latent', z1, 'subspace', w, 'residual', r);
        else
            v = w(:,:,:,:,pg) * z;
        end
        iphi = exponentiateVelocity(v, 'iphi', 'vs', vs, 'prm', prm);
        mu1 = warp(iphi, mu);
        if iscat
            mu1 = reconstructProbaTemplate(mu1);
        end
        mu1 = colorimage(mu1, x, slice);
        image(mu1);
        switch slice
            case 1
                daspect([1/vs(2) 1/vs(3) 1]);
            case 2
                daspect([1/vs(1) 1/vs(3) 1]);
            case 3
                daspect([1/vs(1) 1/vs(2) 1]);
        end
        axis off
        drawnow
        frames(i) = getframe(gcf);
        i = i + 1;
    end
    close(fig);
    
    fprintf('Save frames\n');
    
    if endsWith(fname, 'avi')
        vid = VideoWriter(fname);
        open(vid);
    elseif endsWith(fname, 'mp4')
        vid = VideoWriter(fname, 'MPEG-4');
        open(vid);
    end
    for i=[ceil(n/2):n n:-1:1 1:ceil(n/2)]
        fprintf('z: %f\n', allz(i));
        if endsWith(fname, 'gif')
            [imind, cm] = rgb2ind(frame2im(frames(i)), 256);
            if i == 1
                imwrite(imind, cm, fname, 'gif', 'Loopcount', inf); 
            else
                imwrite(imind, cm, fname, 'gif', 'WriteMode', 'append'); 
            end
        elseif endsWith(fname, 'avi') || endsWith(fname, 'mp4') 
            writeVideo(vid, frames(i));
        end
    end
    if endsWith(fname, 'avi') || endsWith(fname, 'mp4') 
        close(vid);
    end

end

function mu = colorimage(mu, z, slice)
    if nargin < 3
        slice = 3;
    end
    dim = [size(mu) 1 1];
    if nargin < 2 || ~isfinite(z)
        switch slice
            case {2,3}
                z = ceil(dim(slice)/2);
            case 1
                z = ceil(dim(slice)/3);
        end
    end
    switch slice
        case 1
            mu = catToColor(reshape(mu(z,:,:,:), [dim(2) dim(3) dim(4)]));
        case 2
            mu = catToColor(reshape(mu(:,z,:,:), [dim(1) dim(3) dim(4)]));
        case 3
            mu = catToColor(reshape(mu(:,:,z,:), [dim(1) dim(2) dim(4)]));
        otherwise
            error('not handled')
    end
end