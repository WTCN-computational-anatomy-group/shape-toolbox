function [opt, dat, model] = pgra_model_data(opt, dat, model)
% FORMAT [opt, dat, model] = pgra_model_data(opt, dat, model)
%
% Initialise all data handling structures


    % ---------------------------------------------------------------------
    % Compute template dimensions
    % ---------------------------------------------------------------------
    
    % - Default lattice
    if ~isfield(opt.tpl, 'lat')
        if ~isfield(opt.tpl, 'vs')
            opt.tpl.vs = [];
        end
        [opt.tpl.lat, opt.tpl.M, opt.tpl.vs] = autoTemplateDim(dat, opt.tpl.vs);
    end
    if ~isfield(opt.tpl, 'vs')
        if isfield(opt.tpl, 'M')
            opt.tpl.vs = sqrt(sum(opt.tpl.M(1:3,1:3).^2));
        else
            opt.tpl.vs = [1 1 1];
        end
    end
    if ~isfield(opt.tpl, 'M')
        model.tpl.M = eye(4);
        model.tpl.M([1 6 11]) = opt.vs;
        model.tpl.M(1) = -opt.tpl.M(1);
    end
    model.tpl.M = opt.tpl.M;
    
    
    % ---------------------------------------------------------------------
    % Create Data structure
    % ---------------------------------------------------------------------
    
    % I am setting all dimensions to [0 0] and let the first "creator"
    % instance set the appropriate dimensions. I will also add a "cleanup"
    % function at the end which ensures the NIfTI header is correct and
    % only keeps user-defined volumes.
    
    % Model
    % -----
    
    level1 = fieldnames(opt.fnames.model);
    for i=1:numel(level1)
        lv1 = level1{i};
        level2 = fieldnames(opt.fnames.model.(lv1));
        for j=1:numel(level2)
            lv2 = level2{j};
            if ~isfield(model, lv1) || ~isfield(model.(lv1), lv2)
                if opt.ondisk.model.(lv1).(lv2)
                    fname = fullfile(opt.dir.model, opt.fnames.model.(lv1).(lv2));
                    model.(lv1).(lv2) = createFA(fname, [0 0], 'float32');
                else
                    model.(lv1).(lv2) = single([]);
                end
            end
        end
    end
    
    
    % Subjects
    % --------
       
    % Initialise level 1 fields
    level1 = fieldnames(opt.fnames.dat);
    level1 = {level1{:} 'q' 'z'}; % Add q cause it is not in fnames but is needed
    for i=1:numel(level1)
        lv1 = level1{i};
        if ~isfield(dat, lv1)
            [dat.(lv1)] = deal(struct);
        end
    end
    
    % Create background class if needed
    if any(strcmpi(opt.model.name, {'categorical', 'multinomial'}))
        for n=1:numel(dat)
            switch opt.model.nc - size(dat(n).f.f, 4)
                case 0
                    % nothing to do
                case 1
                    % create background class
                        if isempty(opt.dir.dat)
                            fname = {dat(n).f.f.fname};
                            [path,fname,ext] = fileparts(fname{1});
                            fname = fullfile(path, ['bg_' fname ext]);
                        else
                            fname = fullfile(opt.dir.dat, num2str(n), 'background.nii');
                        end
                        createDirectory(fileparts(fname));
                        dim = [size(dat(n).f.f) 1];
                        fa = file_array(fname, dim(1:3), 'float32');
                        nii      = nifti;
                        nii.mat  = dat(n).f.M;
                        nii.mat0 = dat(n).f.M;
                        nii.dat  = fa;
                        create(nii);
                        clear fa
                        bg = ones(size(nii.dat), 'single');
                        for k=1:opt.model.nc-1
                            bg = bg - single(dat(n).f.f(:,:,:,k));
                        end
                        nii.dat(:,:,:) = bg;
                        nii.dat.permission = 'ro';
                        dat(n).f.f = cat(4, dat(n).f.f, nii.dat);
                        clear nii
                otherwise
                    error(['The number of existing classes must either ' ...
                           'be %d or %d'], opt.model.nc, opt.model.nc-1);
            end
        end
    end
    
    % Create only non-existant variables
    level1 = fieldnames(opt.fnames.dat);
    for i=1:numel(level1)
        lv1 = level1{i};
        level2 = fieldnames(opt.fnames.dat.(lv1));
        for j=1:numel(level2)
            lv2 = level2{j};
            for n=1:numel(dat)
                if ~isfield(dat(n).(lv1), lv2)
                    if opt.ondisk.dat.(lv1).(lv2)
                        if isempty(opt.dir.dat)
                            fname = {dat(n).f.f.fname};
                            [path,fname,ext] = fileparts(fname{1});
                            fname = fullfile(path, [opt.fnames.dat.(lv1).(lv2) fname ext]);
                        else
                            fname = fullfile(opt.dir.dat, num2str(n), opt.fnames.dat.(lv1).(lv2));
                        end
                        dat(n).(lv1).(lv2) = createFA(fname, [0 0], 'float32');
                    else
                        dat(n).(lv1).(lv2) = single([]);
                    end
                end
            end
        end
    end
    
    % Create all output directories
    % -----------------------------
    
    level1 = fieldnames(dat);
    for i=1:numel(level1)
        lv1 = level1{i};
        for n=1:numel(dat)
            level2 = fieldnames(dat(n).(lv1));
            for j=1:numel(level2)
                lv2 = level2{j};
                if isa(dat(n).(lv1).(lv2), 'file_array') ...
                        && ~any(strcmpi({dat(n).(lv1).(lv2).permission}, 'ro'))
                    createDirectory(fileparts(dat(n).(lv1).(lv2).fname));
                end
            end
        end
    end
    level1 = fieldnames(model);
    for i=1:numel(level1)
        lv1 = level1{i};
        level2 = fieldnames(model.(lv1));
        for j=1:numel(level2)
            lv2 = level2{j};
            if isa(model.(lv1).(lv2), 'file_array') ...
                    && ~any(strcmpi({model.(lv1).(lv2).permission}, 'ro'))
                createDirectory(fileparts(model.(lv1).(lv2).fname));
            end
        end
    end
    
end

% =========================================================================
% HELPERS
% =========================================================================

function createDirectory(dirpath)
    if ~exist(dirpath, 'dir')
        [parent, dirname] = fileparts(dirpath);
        createDirectory(parent);
        mkdir(parent, dirname);
    end
end

% ===

function f = createFA(fname, dim, dtype, do_create, M)

    if nargin < 5
        M = eye(4);
        if nargin < 4
            do_create = '';
        end
    end

    f      = file_array(fname, dim, dtype);
    n      = nifti;
    n.dat  = f;
    n.mat0 = M;
    f      = n.dat;
    
    if strcmpi(do_create, 'create')
        create(n);
    end
end