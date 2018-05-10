function model = aggregateTemplateGradHess(dat, model, opt)
% FORMAT model = aggregateTemplateGradHess(dat, model, opt)
% dat   - Subject-specific data
% model - Model-specific data
% opt   - Options
%
% All inputs are structures that can either be in memory or on disk in the
% form of a mat file. In the latter case, it is read and, if needed,
% written back.
%--------------------------------------------------------------------------
% Aggregate template gradient and hessian:
%   >> tpl.g / tpl.h
%__________________________________________________________________________
% Copyright (C) 2018 Wellcome Centre for Human Neuroimaging

    % =====================================================================
    % Read input from disk (if needed)
    [dat, ~, model, modelpath, opt] = fileToStruct(dat, model, opt);

    % =====================================================================
    % Prepare on disk
    if isfield(opt, 'fastdir') && ~isempty(opt.fastdir)
        if exist('java.util.UUID', 'class')
            uid = char(java.util.UUID.randomUUID()); % Session UID
        else
            uid = datestr(now, 'yyyymmddTHHMMSS'); % Hopfully no double call in the same second
        end
        mkdir
    	dir = fullfile(opt.fastdir, uid);
        mkdir(dir);
        
        g = file_array(fullfile(dir,'g.nii'));
        h = file_array(fullfile(dir,'h.nii'));
    else
        g = model.tpl.g;
        h = model.tpl.h;
    end
    if iscell(dat)
        dat1 = dat{1};
    else
        dat1 = dat(1);
    end
    if ischar(dat1)
        dat1 = load(dat1);
    end
    g = prepareOnDisk(g,size(dat1.tpl.g));
    h = prepareOnDisk(h,size(dat1.tpl.h));
    for z=1:size(model.tpl.g, 3)
        g(:,:,z,:) = 0;
        h(:,:,z,:) = 0;
    end
    
    % =====================================================================
    % Aggregate data
    for n=1:numel(dat)
        if iscell(dat)
            dat1 = dat{n};
        else
            dat1 = dat(n);
        end
        if ischar(dat1)
            dat1 = load(dat1);
        end
        if defval(dat1.f, '.observed', true)
            if isfield(opt, 'fastdir') && ~isempty(opt.fastdir)
                g1 = rmarray(dat1.tpl.g);
                h1 = rmarray(dat1.tpl.h);
                copyfile(dat1.tpl.g.fname, fullfile(dir,'g1.nii'));
                copyfile(dat1.tpl.h.fname, fullfile(dir,'h1.nii'));
                rmarray(dat1.tpl.g);
                rmarray(dat1.tpl.h);
                g1.fname = fullfile(opt.fastdir,'g1.nii');
                h1.fname = fullfile(opt.fastdir,'h1.nii');
            else
                g1 = dat1.tpl.g;
                h1 = dat1.tpl.h;
            end
            for z=1:size(model.tpl.g, 3)
                g(:,:,z,:) = g(:,:,z,:) + g1(:,:,z,:);
                h(:,:,z,:) = h(:,:,z,:) + h1(:,:,z,:);
            end
            rmarray(g1);
            rmarray(h1);
        end
    end
    
    % =====================================================================
    % Push data
    if isfield(opt, 'fastdir') && ~isempty(opt.fastdir)
        copyfile(g.fname, model.tpl.g.fname);
        copyfile(h.fname, model.tpl.h.fname);
        rmarray(g);
        rmarray(h);
        gname = model.tpl.g.fname;
        hname = model.tpl.h.fname;
        dat.tpl.g = g;
        dat.tpl.g.fname = gname;
        dat.tpl.h = h;
        dat.tpl.h.fname = hname; 
        rmdir(dir, 's');
    else
        model.tpl.g = g;
        model.tpl.h = h;
    end
    
        
    % =====================================================================
    % Exit
    model = structToFile(model, modelpath);
end