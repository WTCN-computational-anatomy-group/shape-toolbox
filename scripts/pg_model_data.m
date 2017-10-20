function [opt, dat, model] = pg_model_data(opt, dat, model)
% FORMAT [opt, dat, model] = pg_model_data(opt, dat, model)
%
% Initialise all data handling structures


   
    % ---------------------------------------------------------------------
    % Model
    % ---------------------------------------------------------------------
    
    % Compute template dimensions
    % ---------------------------
    
    % That bit is quite messy, could be cleaned up.
    
    % - Default lattice
    if ~isfield(opt, 'vs') || ~isfield(opt, 'Mmu') || ~isfield(opt, 'lat')
        vartocheck = {'mu', 'a', 'w'};
        for i=1:numel(vartocheck)
            var = vartocheck{i};
            if isfield(model, var)
                if ~isfield(opt, 'lat')
                    opt.lat = [size(model.(var)) 1 1];
                    opt.lat = opt.lat(1:3);
                end
                if (~isfield(opt, 'vs') || ~isfield(model, 'Mmu')) && ...
                        isa(model.mu, 'file_array') && ...
                        endsWith(model.(var).fname, '.nii')
                    n = nifti(model.(var).fname);
                    model.Mmu = n.mat0;
                    opt.vs = sqrt(sum(model.Mmu(1:3,1:3).^2));
                end
                break
            end
        end
    end
    if ~isfield(opt, 'lat')
        if ~isfield(opt, 'vs')
            opt.vs = [];
        end
        [opt.lat, model.Mmu, opt.vs] = autoTemplateDim(dat, opt.vs);
    end
    if ~isfield(opt, 'vs')
        opt.vs = [1 1 1];
    end
    if ~isfield(model, 'Mmu')
        model.Mmu = eye(4);
        model.Mmu([1 6 11]) = opt.vs;
        model.Mmu(1) = -model.Mmu(1);
    end
    
    
    % Create Data structure
    % ---------------------
    
    % I am supposing that if some data is provided through the model or dat
    % structure, we should respect its 'on/off disk' format, and not try to
    % follow opt.ondisk.
    
    % I am setting all dimensions to [0 0] and let the first "creator"
    % instance set the appropriate dimensions. I will also add a "cleanup"
    % function at the end which ensures the NIfTI header is correct and
    % only keeps user-defined volumes.
    
    % --- Model
    
    modelvar = {'a', 'mu', 'gmu', 'w', 'dw', 'gw', 'hw', ...
                'ww', 'zz', 'z', 'Sz', 'Az'};
            
    for i=1:numel(modelvar)
        var = modelvar{i};
        if ~isfield(model, var)
            if opt.ondisk.model.(var) && ~isempty(opt.fnames.model.(var))
                model.(var) = createFA(opt.fnames.model.(var), [0 0], 'float32');
            else
                model.(var) = single([]);
            end
        end
    end
    if ~isfield(model, 'wpz')
        model.wpz = opt.wpz;
    end
    
    % --- Subjects
            
    datvar = {'wmu', 'iphi', 'ipsi', 'v', 'pf', 'c', 'zz', 'z', 'Sz', ...
              'gv', 'hv', 'g', 'h'};
    
    for j=1:numel(datvar)
        var = datvar{j};
        if ~isfield(dat, var)
            if isfield(opt.ondisk.dat, var) ...
                    && opt.ondisk.dat.(var) ...
                    && isfield(opt.fnames.dat, var) ...
                    && ~isempty(opt.fnames.dat.(var))
                for i=1:opt.N
                    if ~isempty(opt.fnames.dat.(var){i})
                        dat(i).(var) = createFA(opt.fnames.dat.(var){i}, [0 0], 'float32');
                    else
                        dat(i).(var) = single([]);
                    end
                end
            else
                [dat.(var)] = deal(single([]));
            end
        end
    end
    
    if ~isfield(dat, 'llm')
        [dat.llm] = deal(double([]));
    end
    if ~isfield(dat, 'okz')
        [dat.okz] = deal(false);
    end
    
    % Create all output directories
    % -----------------------------
    
    for i=1:opt.N
        fields = fieldnames(dat(i));
        for j=1:numel(fields)
            if isa(dat(i).(fields{j}), 'file_array')
                createDirectory(fileparts(dat(i).(fields{j}).fname));
            end
        end
    end
    fields = fieldnames(model);
    for j=1:numel(fields)
        if isa(model.(fields{j}), 'file_array')
            createDirectory(fileparts(model.(fields{j}).fname));
        end
    end
    
end

% ===

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