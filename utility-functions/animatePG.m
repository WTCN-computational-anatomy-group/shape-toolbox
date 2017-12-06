function frames = animatePG(model, opt, pg, fname)

    if nargin < 4
        fname = '';
    end
        
    allz = -0.3:0.01:0.3;

    Z = ceil(opt.lat(3)/2);
    [X,Y] = ndgrid(1:opt.lat(1), 1:opt.lat(2));
    U = reshape(model.w(:,:,Z,1,5), [opt.lat(1) opt.lat(2)]);
    V = reshape(model.w(:,:,Z,2,5), [opt.lat(1) opt.lat(2)]);
    
    loops = numel(allz);
    frames(loops) = struct('cdata',[],'colormap',[]);
    i = 1;
    for z=allz
        v = model.w(:,:,:,:,pg) * z; % / model.regz(pg,pg);
        iphi = exponentiateVelocity(v, 'iphi', 'vs', opt.vs, 'prm', opt.prm);
        mu = warp(iphi, model.mu);
        mu = colorimage(mu, Z);
        image(mu);
        hold on
        quiver(Y,X,V,U, 'w')
        axis off
        hold off
        drawnow
        frames(i) = getframe(gcf);
        if endsWith(fname, 'gif')
            [imind, cm] = rgb2ind(frame2im(frames(i)), 256);
            if i == 1
                imwrite(imind, cm, fname, 'gif', 'Loopcount', inf); 
            else
                imwrite(imind, cm, fname, 'gif', 'WriteMode', 'append'); 
            end
        elseif endsWith(fname, 'avi')
            if i == 1
                vid = VideoWriter(fname);
                open(vid);
            end
            writeVideo(vid, frames(i));
        end
        i = i + 1;
    end
    
    if endsWith(fname, 'avi')
        for i=numel(frames):-1:1
            writeVideo(vid, frames(i));
        end
        close(vid);
    end
    close

end

function mu = colorimage(mu, z)
    dim = size(mu);
    if nargin < 2
        z = ceil(dim(3)/2);
    end
    mu = reshape(mu(:,:,z,1:3), [dim(1) dim(2) 3]);
end