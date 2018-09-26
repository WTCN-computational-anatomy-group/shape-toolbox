function [opt, dat, model] = shape_input(input, opt)
% FORMAT [opt, dat, model] = shape_input(input, opt)
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
        % Ensure that subjects are along the second dimension
        if size(input.f, 2) == 1
            input.f = input.f';
        end
        % - Build data structure
        f = cell(1, size(input.f, 2));
        M = cell(1, size(input.f, 2));
        for i=1:size(input.f, 2)        % Subjects
            for j=1:size(input.f, 1)    % Classes/Modalities
                if ~exist(input.f{j,i}, 'file')
                    warning('Image %s does not exist. We''ll remove it.', input.f{j,i})
                    f = f(1:end-1);
                    M = M(1:end-1);
                else
                    n = nifti(input.f{j,i});
                    n.dat.permission = 'ro';
                    if numel(size(n.dat)) > 4 && size(n.dat, 4) == 1
                        n.dat.dim = [n.dat.dim(1:3) n.dat.dim(5:end)];
                    end
                    if isa(f{i}, 'file_array')
                        f{i} = cat(4, f{i}, n.dat);
                    else
                        f{i} = n.dat;
                        M{i} = n.mat0;
                    end
                end
            end
        end
    else
        f = {};
        M = {};
    end
    opt.f.N = numel(f);
    if opt.f.N
        if size(f{1}, 3) == 1
            opt.model.dim = 2;
        else
            opt.model.dim = 3;
        end
        if ~isfield(opt, 'model') || ~isfield(opt.model, 'nc')
            opt.model.nc = size(f{1}, 4);
        end
    end
    
    % ---------------------------------------------------------------------
    % Input velocity
    % ---------------------------------------------------------------------

    if isfield(input, 'v') && numel(input.v) > 0
        % - If needed convert from char array to cell array
        if ischar(input.v)
            fnames = cell(1, size(input.v, 1));
            for n=1:numel(fnames)
                fnames{n} = deblank(input.v(n,:));
            end
            input.v = fnames;
        end
        % - Build data structure
        v = cell(1, numel(input.v));
        for i=1:numel(input.v)
            if ~exist(input.v{i}, 'file')
                warning('Velocity %s does not exist. We''ll remove it.', input.v{i})
                v = v(1:end-1);
            else
                n = nifti(input.v{i});
                n.dat.permission = 'ro';
                v{i} = n.dat;
                if numel(size(v{i})) > 4 && size(v{i}, 4) == 1
                    v{i}.dim = [v{i}.dim(1:3) v{i}.dim(5:end)];
                end
                if ~isfield(opt, 'tpl') || ~isfield(opt.tpl, 'lat')
                    opt.tpl.lat = v{i}.dim(1:3);
                    opt.tpl.M   = n.mat0;
                    opt.tpl.vs  = sqrt(sum(n.mat0(1:3,1:3).^2));
                end
            end
        end
    else
        v = {};
    end
    opt.v.N = numel(v);
    if opt.v.N
        if size(v{1}, 3) == 1
            opt.model.dim = 2;
        else
            opt.model.dim = 3;
        end
    end
    
    % ---------------------------------------------------------------------
    % Build data structure
    % ---------------------------------------------------------------------
    
    dat = struct;
    % images
    [dat(1:opt.f.N+opt.v.N).f] = deal(struct);
    for n=1:opt.f.N
        dat(n).f = struct('observed', true, 'f', f{n}, 'M', M{n});
    end
    [dat(opt.f.N+1:end).f]     = deal(struct('observed', false));
    % velocities
    [dat(1:opt.f.N+opt.v.N).v] = deal(struct);
    i = 0;
    for n=opt.f.N+1:numel(dat)
        i = i + 1;
        dat(n).v = struct('observed', true, 'v', v{i});
    end
    [dat(1:opt.f.N).v] = deal(struct('observed', false));
    
    
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
