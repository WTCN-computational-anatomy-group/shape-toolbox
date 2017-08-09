function a = setNifti(value, varargin)
% FORMAT a = setNifti(value, (a), 'OptKey', opt_value, ...)
% FORMAT a = setNifti(value, (a), OptStruct)
%
% value  - New value (array, nifti, path)
% (a)    - Existing nifti object
% OptKey -  * 'replace'
%           * 'fname'
%           * 'dtype'
%           * 'dim'
%           * 'value'
% OptStruct - Structure with fields in OptKey

    % --- Create new object if none passed
    if ischar(varargin{1}) || isstruct(varargin{1})
        fa = file_array;
        fa.dtype = 'float32';
        fa.fname = 'foo.nii';
        a = nifti;
        a.dat = fa;
    else
        a = varargin{1};
        varargin = varargin(2:end);
    end
    
    % --- First parsing to check if the value is passed with key 'value'
    if isnumeric(value) && isempty(value) && ~isempty(varargin)
        if ischar(varargin{1})
            for i=1:numel(varargin)
                if ischar(varargin{i}) && strcmpi(varargin{i}, 'value')
                    value = varargin{i+1};
                    varargin = [varargin(1:i-1) varargin(1+3:end)];
                    break
                end
            end
        elseif isstruct(varargin{1})
            if isfield(varargin{1}, 'value')
                value = varargin{1}.value;
                varargin{1} = rmfield(varargin{1}, 'value');
            end
        end
    end

    
    % --- If path provided, read it
    if ischar(value)
        value = nifti(value);
    end
    
    % --- Set default arguments
    opt = struct('replace', true);
    if isa(value, 'nifti')
        opt.fname   = value.dat.fname;
        opt.dim     = value.dat.dim;
        opt.dtype   = value.dat.dtype;
    elseif ~isempty(value)
        opt.dim     = size(value);
        opt.dtype   = matlab2FileArrayType(class(value));
    end 
    
    % --- Read optional arguments
    if ~isempty(varargin) && isstruct(varargin{1})
        opt = update_options(opt, varargin{1});
    else
        opt = parse_varargin(varargin, opt);
    end
    
    % --- Eventually, completely replace the nifti object
    if opt.replace && isa(value, 'nifti')
        a = value;
    end
    
    % --- Apply optional changes
    if isfield(opt, 'fname')
        a.dat.fname = opt.fname;
    end
    if isfield(opt, 'dim')
        a.dat.dim = opt.dim;
    end
    if isfield(opt, 'dtype')
        a.dat.dtype = opt.dtype;
    end
    if isfield(opt, 'directory')
        [~, name, ext] = fileparts(a.dat.fname);
        a.dat.fname = [opt.directory filesep name ext];
    end
    
    % --- Transfer data
    if ~(opt.replace && isa(value, 'nifti')) && ...
       ~(isnumeric(value) && isempty(value))
        if ~exist(a.dat.fname, 'file')
            create(a);
        end
        a.dat(:) = value(:);
    end
    
    % --- Convert relative path to absolute path
    if isfield(a.dat, 'fname')
        if isunix && ~startsWith(a.fname, '/')
            a.fname = [pwd a.dat.fname];
        elseif ispc && isempty(regexpi(a.fname, '^[A-Z]+:[\\/]', 'start'))
            a.fname = [pwd a.dat.fname];
        end
    end
end