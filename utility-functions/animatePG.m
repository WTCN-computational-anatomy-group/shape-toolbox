function frames = animatePG(model, opt, pg, fname, sigma)

    if nargin < 5
        sigma = 3;
        if nargin < 4
            fname = '';
        end
    end
    
    cov = inv(model.Az);
    nf = 21;
    allz = linspace(-sigma*cov(pg,pg), sigma*cov(pg,pg), nf);

    Z = ceil(opt.lat(3)/2);
    
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
        mu = colorimage(mu, Z);
        image(mu);
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

function mu = colorimage(mu, z)
    dim = [size(mu) 1 1];
    if nargin < 2
        z = ceil(dim(3)/2);
    end
    mu = catToColor(reshape(mu(:,:,z,:), [dim(1) dim(2) dim(4)]));
end