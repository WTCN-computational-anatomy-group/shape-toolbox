function pgra_plot_all(model, opt)
% Plot PG + lower bound stuff
% This function is highly specific to this particular model. Not sure I can
% come up with a generic "plotModel" function, even though it could be
% nice.
    
    if opt.ui.verbose
        if islogical(opt.ui.ftrack) && ~opt.ui.ftrack
            return
        end
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
        if dim(4) > 1
            tpl = reshape(tpl, [dim(1:2) dim(4)]);
        end
        tpl = permute(tpl, [2 1 3]);
        asp = 1./[opt.tpl.vs(2) opt.tpl.vs(1) 1];
        image(tpl(end:-1:1,:,:));
        daspect(asp);
        axis off
        npg = npg - 1;
        if size(model.tpl.mu, 3) > 1
            title('template (axial)')
            i = i + 1;
            subplot(nh, nw, i)
            tpl = model.tpl.mu(:,ceil(size(model.tpl.mu,2)/2),:,:);
            tpl = reshape(tpl, [size(tpl,1) size(tpl,3) size(tpl,4)]);
            tpl = catToColor(tpl);
            tpl = permute(tpl, [2 1 3]);
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
        if isfield(model.lb, 'm')
            subplot(nh,nw,i)
            plot([model.lb.lb.it model.lb.lb.curit], model.lb.m.list, ...
                 colors(mod(i, length(colors))+1))
            title(model.lb.m.name)
        end
        % PG prior
        i = i + 1;
        if isfield(model.lb, 'w')
            subplot(nh,nw,i)
            plot([model.lb.lb.it model.lb.lb.curit], model.lb.w.list, ...
                 colors(mod(i, length(colors))+1))
            title(model.lb.w.name)
        end
        
        % --- Line 3
        
        i = i + 1;
        if opt.q.Mr
            % KL affine
            if isfield(model.lb, 'q')
                subplot(nh,nw,i)
                plot([model.lb.lb.it model.lb.lb.curit], model.lb.q.list, ...
                     colors(mod(i, length(colors))+1))
                title(model.lb.q.name)
            end
        else
            % LL geodesic
            if isfield(model.lb, 'g')
                subplot(nh,nw,i)
                plot([model.lb.lb.it model.lb.lb.curit], model.lb.g.list, ...
                     colors(mod(i, length(colors))+1))
                title(model.lb.g.name)
            end
        end
        % KL residual
        i = i + 1;
        if isfield(model.lb, 'r')
            subplot(nh,nw,i)
            plot([model.lb.lb.it model.lb.lb.curit], model.lb.r.list, ...
                 colors(mod(i, length(colors))+1))
            title(model.lb.r.name)
        end
        % KL latent
        i = i + 1;
        if isfield(model.lb, 'z')
            subplot(nh,nw,i)
            plot([model.lb.lb.it model.lb.lb.curit], model.lb.z.list, ...
                 colors(mod(i, length(colors))+1))
            title(model.lb.z.name)
        end
        
        
        % --- Line 4
        
        i = i + 1;
        if opt.q.Mr
            % KL affine precision
            if isfield(model.lb, 'Aq')
                subplot(nh,nw,i)
                plot([model.lb.lb.it model.lb.lb.curit], model.lb.Aq.list, ...
                     colors(mod(i, length(colors))+1))
                title(model.lb.Aq.name)
            end
        end
        % KL residual precision
        i = i + 1;
        if isfield(model.lb, 'l')
            subplot(nh,nw,i)
            plot([model.lb.lb.it model.lb.lb.curit], model.lb.l.list, ...
                 colors(mod(i, length(colors))+1))
            title(model.lb.l.name)
        end
        % KL latent precision
        i = i + 1;
        if isfield(model.lb, 'Az')
            subplot(nh,nw,i)
            plot([model.lb.lb.it model.lb.lb.curit], model.lb.Az.list, ...
                 colors(mod(i, length(colors))+1))
            title(model.lb.Az.name)
        end
        
        % --- Line 5
        
        % WW
        i = i + 1;
        if isfield(model.pg, 'ww')
            subplot(nh,nw,i)
            imagesc(model.pg.ww), colorbar;
            title('W''LW')
        end
        % Sample covariance
        i = i + 1;
        if isfield(model.z, 'zz')
            subplot(nh,nw,i)
            imagesc(model.z.S + model.z.zz), colorbar;
            title('Sample covariance E[ZZ]')
        end
        % Precision matrix
        i = i + 1;
        if isfield(model.z, 'A')
            subplot(nh,nw,i)
            imagesc(model.z.A), colorbar;
            title('Latent precision E[A]')
        end
        
        drawnow
    end

end