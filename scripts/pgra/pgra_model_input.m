function [opt, dat, model] = pgra_model_input(input, opt)
% FORMAT [opt, dat, model] = pgra_model_input(input, opt)
%
% Parse provided observed images and start setting the data structure.
% Should be called first.

    % ---------------------------------------------------------------------
    % Input images
    % ---------------------------------------------------------------------
    
    if isfield(input, 'f') && numel(input.f) > 0
        % - If needed convert from char array to cell array
        if ischar(input.f)
            fnames = cell(1, size(input.f, 1));
            for n=1:numel(fnames)
                fnames{n} = deblank(input.f(n,:));
            end
            input.f = fnames;
        end
        % - Build data structure
        f = cell(1, numel(input.f));
        M = cell(1, numel(input.f));
        for i=1:numel(input.f)
            if ~exist(input.f{i}, 'file')
                warning('Image %s does not exist. We''ll remove it.', input.f{i})
                f = f(1:end-1);
                M = M(1:end-1);
            else
                n = nifti(input.f{i});
                n.dat.permission = 'ro';
                f{i} = n.dat;
                M{i} = n.mat0;
                if numel(size(f{i})) > 4 && size(f{i}, 4) == 1
                    f{i}.dim = [f{i}.dim(1:3) f{i}.dim(5:end)];
                end
            end
        end
    else
        f = {};
        M = {};
    end
    opt.N = numel(f);
    if opt.N
        if size(f{1}, 3) == 1
            opt.model.dim = 2;
        else
            opt.model.dim = 3;
        end
        opt.model.nc = size(f{1}, 4);
    end
    
    % ---------------------------------------------------------------------
    % Build data structure
    % ---------------------------------------------------------------------
    
    dat = struct;
    [dat(1:opt.N).f] = deal(struct);
    for n=1:opt.N
        dat(n).f = struct('observed', true, 'f', f{n}, 'M', M{n});
    end
    
    % ---------------------------------------------------------------------
    % Build model structure
    % ---------------------------------------------------------------------
    
    model     = struct;
    model.pg  = struct;
    model.tpl = struct;
    
    opt.pg.provided = false;
    if isfield(input, 'w') && numel(input.w) > 0
        if ~exist(input.w, 'file')
            warning('Subspace %s does not exist. We''ll remove it.', input.w)
        else
            n = nifti(input.w);
            n.dat.permission = 'ro';
            w = n.dat;
            if numel(size(w)) > 4 && size(w, 4) == 1
                w.dim = [w.dim(1:3) w.dim(5:end)];
            end
            model.pg.w = w;
            opt.pg.provided = true;
            opt.tpl.lat = w.dim(1:3);
            opt.tpl.M   = n.mat0;
            opt.tpl.vs  = sqrt(sum(n.mat0(1:3,1:3).^2));
        end
    end
    
    opt.tpl.provided = false;
    if isfield(input, 'a') && numel(input.a) > 0
        if ~exist(input.a, 'file')
            warning('Log-template %s does not exist. We''ll remove it.', input.a)
        else
            n = nifti(input.a);
            n.dat.permission = 'ro';
            a = n.dat;
            if numel(size(a)) > 4 && size(a, 4) == 1
                a.dim = [a.dim(1:3) a.dim(5:end)];
            end
            model.tpl.a = a;
            opt.tpl.provided = true;
            opt.tpl.lat = a.dim(1:3);
            opt.tpl.M   = n.mat0;
            opt.tpl.vs  = sqrt(sum(n.mat0(1:3,1:3).^2));
        end
    end
    
    opt.tpl.provided = false;
    if isfield(input, 'mu') && numel(input.mu) > 0
        if ~exist(input.mu, 'file')
            warning('Template %s does not exist. We''ll remove it.', input.mu)
        else
            n = nifti(input.mu);
            n.dat.permission = 'ro';
            mu = n.dat;
            if numel(size(mu)) > 4 && size(mu, 4) == 1
                mu.dim = [mu.dim(1:3) mu.dim(5:end)];
            end
            model.tpl.mu = mu;
            opt.tpl.provided = true;
            opt.tpl.lat = mu.dim(1:3);
            opt.tpl.M   = n.mat0;
            opt.tpl.vs  = sqrt(sum(n.mat0(1:3,1:3).^2));
        end
    end
    
end