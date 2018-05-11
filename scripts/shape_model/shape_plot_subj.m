function ind = shape_plot_subj(dat, model, opt, ind)
% Plot PG + lower bound stuff
% This function is highly specific to this particular model. Not sure I can
% come up with a generic "plotModel" function, even though it could be
% nice.
    
    if opt.ui.verbose && ischar(opt.ui.figure_sub)
        
        if nargin < 4 || isempty(ind)
            N = min(opt.f.N, 5);
            ind = randperm(numel(dat), N);
        else
            N = numel(ind);
        end
        
        f = findobj('Type', 'Figure', 'Name', opt.ui.figure_sub);
        if isempty(f)
            f = figure('Name', opt.ui.figure_sub, 'NumberTitle', 'off');
        end
        set(0, 'CurrentFigure', f);   
        clf(f);
        
        nw = 4;
        nh = N;
        i  = 0;
        
        for j=1:numel(ind)
            i = plot_one(dat(ind(j)), model, opt, i, nh, nw);
        end
        
    end
        
end
    
function i = plot_one(dat, model, opt, i, nh, nw)

    if dat.f.observed

        % -------------------------------------------------------------
        % Observed image
        i = i+1;
        subplot(nh, nw, i);
        im = catToColor(dat.f.f(:,:,ceil(size(dat.f.f,3)/2),:));
        dim = [size(im) 1 1];
        if dim(4) > 1
            im = reshape(im, [dim(1:2) dim(4)]);
        end
        im = permute(im, [2 1 3]);
        vs = sqrt(sum(dat.f.M(1:3,1:3).^2));
        asp = 1./[vs(2) vs(1) 1];
        image(im(end:-1:1,:,:));
        daspect(asp);
        axis off

        % -------------------------------------------------------------
        % Warped template
        i = i+1;
        subplot(nh, nw, i);
        im = catToColor(dat.tpl.wmu(:,:,ceil(size(dat.tpl.wmu,3)/2),:));
        dim = [size(im) 1 1];
        if dim(4) > 1
            im = reshape(im, [dim(1:2) dim(4)]);
        end
        im = permute(im, [2 1 3]);
        vs = sqrt(sum(dat.f.M(1:3,1:3).^2));
        asp = 1./[vs(2) vs(1) 1];
        image(im(end:-1:1,:,:));
        daspect(asp);
        axis off

        % -------------------------------------------------------------
        % Deformation
        i = i + 1;
        subplot(nh,nw,i)
        def = numeric(dat.v.ipsi);
        id = reconstructIPsi(dat.q.A, spm_warps('identity', opt.tpl.lat));
        def = def - id;
        clear id
        im = defToColor(def(:,:,ceil(size(def,3)/2),:));
        clear def
        dim = [size(im) 1 1];
        im = permute(reshape(im, [dim(1:2) dim(4)]), [2 1 3]);
        vs = sqrt(sum(dat.f.M(1:3,1:3).^2));
        asp = 1./[vs(2) vs(1) 1];
        image(im(end:-1:1,:,:));
        daspect(asp);
        axis off

        % -------------------------------------------------------------
        % Latent coordinates
        i = i + 1;
        subplot(nh,nw,i)
        scatter(model.z.Z(1,:),model.z.Z(2,:), 'filled')
        hold on
        scatter(dat.z.z(1), dat.z.z(2), 500, 'red', 'p', 'filled')
        hold off

    else
    end
        
end