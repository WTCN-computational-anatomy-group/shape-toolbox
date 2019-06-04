function shape_plot_all(model, opt)
% Plot PG + lower bound stuff
% This function is highly specific to this particular model. Not sure I can
% come up with a generic "plotModel" function, even though it could be
% nice.
    
    if opt.ui.verbose && ischar(opt.ui.figure_pop)
        f = findobj('Type', 'Figure', 'Name', opt.ui.figure_pop);
        if isempty(f)
            f = figure('Name', opt.ui.figure_pop, 'NumberTitle', 'off');
        end
        set(0, 'CurrentFigure', f);   
        clf(f);
        
        nw = 3;
        nh = 5;
        i  = 0;
        colors = ['b', 'g', 'r', 'c', 'm', 'k'];
        npg = nw;
        
        % --- Line 1
        
        % Template
        if opt.f.N && opt.optimise.tpl.a
            i = i + 1;
            subplot(nh, nw, i)
            if any(strcmpi(opt.model.name, {'categorical', 'multinomial', 'cat', 'mult'}))
                tpl = catToColor(model.tpl.mu(:,:,ceil(size(model.tpl.mu,3)/2),:));
                implot = @(X) image(X);
            else
                tpl = model.tpl.mu(:,:,ceil(size(model.tpl.mu,3)/2),1);
                implot = @(X) colormap(get(imagesc(X), 'Parent'), 'gray');
            end
            dim = [size(tpl) 1 1];
            if dim(4) > 1
                tpl = reshape(tpl, [dim(1:2) dim(4)]);
            end
            tpl = permute(tpl, [2 1 3]);
            asp = 1./[opt.tpl.vs(2) opt.tpl.vs(1) 1];
            implot(tpl(end:-1:1,:,:));
            daspect(asp);
            axis off
            npg = npg - 1;
            if size(model.tpl.mu, 3) > 1
                title('template (axial)')
                i = i + 1;
                subplot(nh, nw, i)
                if any(strcmpi(opt.model.name, {'categorical', 'multinomial', 'cat', 'mult'}))
                    tpl = catToColor(model.tpl.mu(:,ceil(size(model.tpl.mu,2)/2),:,:));
                    implot = @(X) image(X);
                else
                    tpl = model.tpl.mu(:,:,ceil(size(model.tpl.mu,3)/2),1);
                    implot = @(X) colormap(get(imagesc(X), 'Parent'), 'gray');
                end
                dim = [size(tpl) 1 1];
                if dim(4) > 1
                    tpl = reshape(tpl, [dim(1) dim(3) dim(4)]);
                end
                tpl = permute(tpl, [2 1 3]);
                asp = 1./[opt.tpl.vs(3) opt.tpl.vs(1) 1];
                implot(tpl(end:-1:1,:,:));
                daspect(asp);
                axis off
                title('template (coronal)')
                npg = npg - 1;
            else
                title('template')
            end
        end 
        if opt.optimise.pg.w
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
        else
            i = i + npg;
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
            plot([model.lb.lb.it model.lb.lb.curit], model.lb.f.list, ...
                 colors(mod(i, length(colors))+1))
            title(model.lb.f.name)
        else
            i = i +1;
        end
        if opt.optimise.pg.w
            % PG prior
            i = i + 1;
            subplot(nh,nw,i)
            plot([model.lb.lb.it model.lb.lb.curit], model.lb.w.list, ...
                 colors(mod(i, length(colors))+1))
            title(model.lb.w.name)
        else
            i = i +1;
        end
        
        % --- Line 3
        
        if opt.q.Mr && opt.optimise.q.q
            % KL affine
            i = i + 1;
            subplot(nh,nw,i)
            plot([model.lb.lb.it model.lb.lb.curit], model.lb.q.list, ...
                 colors(mod(i, length(colors))+1))
            title(model.lb.q.name)
        elseif opt.v.N && opt.optimise.v.v
            % LL velocity
            i = i + 1;
            subplot(nh,nw,i)
            plot([model.lb.lb.it model.lb.lb.curit], model.lb.v1.list, ...
                 colors(mod(i, length(colors))+1))
            title(model.lb.v1.name)
        elseif opt.optimise.tpl.a && strcmpi(opt.tpl.update, 'map') && sum(opt.tpl.prm) > 0
            % Template prior
            i = i + 1;
            subplot(nh,nw,i)
            plot([model.lb.lb.it model.lb.lb.curit], model.lb.a.list, ...
                 colors(mod(i, length(colors))+1))
            title(model.lb.a.name)
        else
            i = i + 1;
        end
        if opt.f.N && opt.optimise.v.v
            % KL velocity
            i = i + 1;
            subplot(nh,nw,i)
            plot([model.lb.lb.it model.lb.lb.curit], model.lb.v2.list, ...
                 colors(mod(i, length(colors))+1))
            title(model.lb.v2.name)
        else
            i = i + 1;
        end
        if opt.optimise.z.z
            % KL latent
            i = i + 1;
            subplot(nh,nw,i)
            plot([model.lb.lb.it model.lb.lb.curit], model.lb.z.list, ...
                 colors(mod(i, length(colors))+1))
            title(model.lb.z.name)
        else
            i = i + 1;
        end
        
        
        % --- Line 4
        
        if opt.q.Mr && opt.optimise.q.A
            % KL affine precision
            i = i + 1;
            subplot(nh,nw,i)
            plot([model.lb.lb.it model.lb.lb.curit], model.lb.Aq.list, ...
                 colors(mod(i, length(colors))+1))
            title(model.lb.Aq.name)
        else
            i = i + 1;
        end
        if opt.optimise.v.l
            % KL residual precision
            i = i + 1;
            subplot(nh,nw,i)
            plot([model.lb.lb.it model.lb.lb.curit], model.lb.l.list, ...
                 colors(mod(i, length(colors))+1))
            title(model.lb.l.name)
        else
            i = i + 1;
        end
        if opt.optimise.z.A
            % KL latent precision
            i = i + 1;
            subplot(nh,nw,i)
            plot([model.lb.lb.it model.lb.lb.curit], model.lb.Az.list, ...
                 colors(mod(i, length(colors))+1))
            title(model.lb.Az.name)
        else
            i = i + 1;
        end
        
        % --- Line 5
        
        if opt.optimise.pg.w
            % WW
            i = i + 1;
            subplot(nh,nw,i)
            imagesc(model.pg.ww), colorbar;
            title('W''LW')
        else
            i = i + 1;
        end
        if opt.optimise.z.z
            % Sample covariance
            i = i + 1;
            subplot(nh,nw,i)
            imagesc(model.z.S + model.z.zz), colorbar;
            title('Sample covariance E[ZZ]')
        else
            i = i + 1;
        end
        if opt.optimise.z.A || opt.optimise.z.z
            % Precision matrix
            i = i + 1;
            subplot(nh,nw,i)
            imagesc(model.z.A), colorbar;
            title('Latent precision E[A]')
        else
            i = i + 1;
        end
        
        drawnow
    end

end