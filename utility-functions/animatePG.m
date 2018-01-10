function frames = animatePG(model, opt, pg, fname, nsigma, dir, x)
% FORMAT (frames) = animatePG(model, opt, pg, (fname), (nsigma))
% model  - Model obtained from pgra_model
% opt    - Options obtained from pgra_model
% pg     - Principal geodesic along which to shoot (1, 2, ..., K)
% fname  - Filename of the output video (.avi, .gif) [do not write video]
% nsigma - Shooting range in terms of number of standard deviation [3]
% frames - Structure array of Matlab frames. 
%          They can be given to the movie function.

    if nargin < 7
        x = nan;
        if nargin < 6
            dir = 3;
            if nargin < 5
                nsigma = 3;
                if nargin < 4
                    fname = '';
                end
            end
        end
    end
    
    cov  = inv(model.Az);
    sd   = sqrt(cov(pg,pg));
    nf   = 11; % Number of frames to compute (should be a parameter really)
    allz = linspace(-nsigma*sd, nsigma*sd, nf);
    
    loops = numel(allz);
    frames(loops) = struct('cdata',[],'colormap',[]);
    i = 1;
    fprintf('Build and record frames\n');
    fig = figure;
    for z=allz
        fprintf('z: %f\n', z);
        v = model.w(:,:,:,:,pg) * z;
        iphi = exponentiateVelocity(v, 'iphi', 'vs', opt.vs, 'prm', opt.prm);
        mu = warp(iphi, model.mu);
        mu = colorimage(mu, x, dir);
        image(mu);
        switch dir
            case 1
                daspect([1/opt.vs(2) 1/opt.vs(3) 1]);
            case 2
                daspect([1/opt.vs(1) 1/opt.vs(3) 1]);
            case 3
                daspect([1/opt.vs(1) 1/opt.vs(2) 1]);
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
    end
    for i=[ceil(nf/2):nf nf:-1:1 1:ceil(nf/2)]
        fprintf('z: %f\n', allz(i));
        if endsWith(fname, 'gif')
            [imind, cm] = rgb2ind(frame2im(frames(i)), 256);
            if i == 1
                imwrite(imind, cm, fname, 'gif', 'Loopcount', inf); 
            else
                imwrite(imind, cm, fname, 'gif', 'WriteMode', 'append'); 
            end
        elseif endsWith(fname, 'avi')
            writeVideo(vid, frames(i));
        end
    end
    if endsWith(fname, 'avi')
        close(vid);
    end

end

function mu = colorimage(mu, z, dir)
    if nargin < 3
        dir = 3;
    end
    dim = [size(mu) 1 1];
    if nargin < 2 || ~isfinite(z)
        switch dir
            case {2,3}
                z = ceil(dim(dir)/2);
            case 1
                z = ceil(dim(dir)/3);
        end
    end
    switch dir
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