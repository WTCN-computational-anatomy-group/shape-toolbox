function mu = updateMuCategoricalML(varargin)
%__________________________________________________________________________
%
% Closed form M-step update of the template for the Categorical model.
% (Maximum likelihood update: no prior on Mu)
%
% -------------------------------------------------------------------------
%
% FORMAT mu = updateMuCategoricalML('f', ..., 'c', ..., ('bb', ...), ...)
%
% REQUIRED KEYWORD (LIST)
% -----------------------
% f   - Images pushed in template space
% c   - Pushed voxel counts
%
% OPTIONAL KEYWORD (LIST)
% -----------------------
% bb   - Bounding boxes (if different from template)
%
% KEYWORD ARGUMENTS
% -----------------
% lat  - Template lattice [first image]
% fwhm - Smoothing kernel used as pseudo prior [do not use]
% loop - How to split processing: 'component', 'slice', 'none' or '' [auto]
% par  - Distribute compute [auto]
%
% OUTPUT
% ------
% mu   - Updated template
%__________________________________________________________________________
    
    i = 1;
    if ischar(varargin{1})
        i = i + 1;
        N = 0;
        while i <= numel(varargin) && ~ischar(varargin{i})
            i = i + 1;
            N = N + 1;
        end
    end
    
    f = {};
    if ischar(varargin{1}) && strcmpi(varargin{1}, 'f')
        f = varargin(2:N+1);
        varargin = varargin(N+2:end);
    end
    c = {};
    if ischar(varargin{1}) && strcmpi(varargin{1}, 'c')
        c = varargin(2:N+1);
        varargin = varargin(N+2:end);
    end
    bb = {};
    if ischar(varargin{1}) && strcmpi(varargin{1}, 'bb')
        bb = varargin(2:N+1);
        varargin = varargin(N+2:end);
    end

    if numel(c) ~= N
        error('There should be as many count as intensity images')
    end
    if numel(bb) > 0 && numel(bb) ~= N
        error('There should be as many bounding boxes as intensity images')
    end

    % --- Parse inputs
    p = inputParser;
    p.FunctionName = 'updateMuCategoricalML';
    p.addParameter('lat',    [],    @isnumeric);
    p.addParameter('fwhm',   0,     @isnumeric);
    p.addParameter('loop',   '',    @(X) ischar(X) && any(strcmpi(X, {'slice', 'subject', 'none', ''})));
    p.addParameter('par',    false, @isscalar);
    p.addParameter('debug',  false, @isscalar);
    p.addParameter('output', false);
    p.parse(varargin{:});
    lat    = p.Results.lat;
    fwhm   = p.Results.fwhm;
    loop   = p.Results.loop;
    par    = p.Results.par;
    debug  = p.Results.debug;
    output = p.Results.output;
    
    if debug, fprintf('* updateMuCategoricalML\n'); end

    if isempty(lat)
        error('For now, the lattice size MUST be provided')
    end

    [par, loop] = autoParLoop(par, loop, isa(f{1}, 'file_array'), size(f{1}, 3));
    if fwhm > 0 && strcmpi(loop, 'slice')
        loop = 'none';
    end
    
    switch lower(loop)
        case 'none'
            if debug, fprintf('   - No loop\n'); end
            mu = loopNone(f, c, bb, lat, output, fwhm);
        case 'slice'
            if debug
                if par > 0
                    fprintf('   - Parallelise on slices\n');
                else
                    fprintf('   - Serialise on slices\n');
                end
            end
            mu = loopSlice(f, c, bb, lat, par, output);
        otherwise
            error('Unknown loop type ''%s''', loop)
    end

end

function mu = loopSlice(f, c, bb, lat, par, output)

    nc = size(f{1}, 4);
    mu = prepareOnDisk(output, [lat nc], 'type', 'float32');
    
    % --- Compute count
    tmpc = zeros(lat, 'single');
    for n=1:numel(f)
        if ~isempty(bb)
            bx = bb{n}.x;
            by = bb{n}.y;
            bz = bb{n}.z;
        else
            bx = 1:lat(1);
            by = 1:lat(2);
            bz = 1:lat(3);
        end
        tmpc(bx,by,bz) = tmpc(bx,by,bz) + single(numeric(c{n}));
    end
    
    % --- Compute mu
    if ~par
        for z=1:lat(3)
            tmpf = zeros([lat(1:2) 1 nc], 'single');
            for n=1:numel(f)
                if ~isempty(bb)
                    bx = bb{n}.x;
                    by = bb{n}.y;
                    bz = bb{n}.z;
                else
                    bx = 1:lat(1);
                    by = 1:lat(2);
                    bz = 1:lat(3);
                end
                fz = z-bz(1)+1;
                if fz >= 1 && fz <= size(f{n}, 3)
                    tmpf(bx,by,1,:) = tmpf(bx,by,1,:) + single(f{n}(:,:,fz,:));
                end
            end
            tmpf = bsxfun(@rdivide, tmpf, tmpc(:,:,z));
            tmpf = max(tmpf, eps('single'));
            tmpf = bsxfun(@rdivide, tmpf, tmpf(:,:,:,nc));
            mu(:,:,z,:) = log(tmpf);
        end
    elseif isa(f{1}, 'file_array')
        parfor (z=1:lat(3), par)
            tmpf = zeros([lat(1:2) 1 nc], 'single');
            for n=1:numel(f)
                if ~isempty(bb)
                    bx = bb{n}.x;
                    by = bb{n}.y;
                    bz = bb{n}.z;
                else
                    bx = 1:lat(1);
                    by = 1:lat(2);
                    bz = 1:lat(3);
                end
                fz = z-bz(1)+1;
                if fz >= 1 && fz <= size(f{n}, 3)
                    tmpf(bx,by,1,:) = tmpf(bx,by,1,:) + single(slicevol(f{n}, fz, 3));
                end
            end
            tmpf = bsxfun(@rdivide, tmpf, tmpc(:,:,z));
            tmpf = max(tmpf, eps('single'));
            tmpf = bsxfun(@rdivide, tmpf, tmpf(:,:,:,nc));
            mu(:,:,z,:) = log(tmpf);
        end
    else
        parfor (z=1:lat(3), par)
            tmpf = zeros([lat(1:2) 1 nc], 'single');
            for n=1:numel(f)
                if ~isempty(bb)
                    bx = bb{n}.x;
                    by = bb{n}.y;
                    bz = bb{n}.z;
                else
                    bx = 1:lat(1);
                    by = 1:lat(2);
                    bz = 1:lat(3);
                end
                fz = z-bz(1)+1;
                if fz >= 1 && fz <= size(f{n}, 3)
                    tmpf(bx,by,1,:) = tmpf(bx,by,1,:) + single(f{n}(:,:,fz,:));
                end
            end
            tmpf = bsxfun(@rdivide, tmpf, tmpc(:,:,z));
            tmpf = max(tmpf, eps('single'));
            tmpf = bsxfun(@rdivide, tmpf, tmpf(:,:,:,nc));
            mu(:,:,z,:) = log(tmpf);
        end
    end
    
    if ~isempty(output)
        mu = saveOnDisk(output, mu, 'name', 'mu');
    end
end

function mu = loopNone(f, c, bb, lat, output, fwhm)

    nc = size(f{1}, 4);
    
    mu   = zeros([lat nc], 'single');
    tmpc = zeros(lat, 'single');
    for n=1:numel(f)
        if ~isempty(bb)
            bx = bb{n}.x;
            by = bb{n}.y;
            bz = bb{n}.z;
        else
            bx = 1:lat(1);
            by = 1:lat(2);
            bz = 1:lat(3);
        end
        mu(bx,by,bz,:) = mu(bx,by,bz,:) + single(numeric(f{n}));
        tmpc(bx,by,bz) = tmpc(bx,by,bz) + single(numeric(c{n}));
    end
    mu   = smooth_gaussian(mu, fwhm);
    tmpc = smooth_gaussian(tmpc, fwhm);
    mu   = bsxfun(@rdivide, mu, tmpc);
    mu   = max(mu, eps('single'));
    mu   = bsxfun(@rdivide, mu, mu(:,:,:,nc));
    mu   = log(mu);
    
    if ~isempty(output)
        mu = saveOnDisk(output, mu, 'name', 'mu');
    end
end