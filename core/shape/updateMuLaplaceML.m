function mu = updateMuLaplaceML(varargin)
%__________________________________________________________________________
%
% Closed form M-step update of the template for the Laplace noise model.
% (Maximum likelihood update: no prior on Mu)
%
% -------------------------------------------------------------------------
%
% FORMAT mu = updateMuNormalML('f', ..., 'c', ...,
%                              ('b', ...),  ('bb', ...), ...)
%
% REQUIRED KEYWORD (LIST)
% -----------------------
% f   - Images pushed in template space
% c   - Pushed voxel counts
%
% OPTIONAL KEYWORD (LIST)
% -----------------------
% b    - Individual noise variances
% bb   - Bounding boxes (if different from template)
%
% KEYWORD ARGUMENTS
% -----------------
% lat  - Template lattice [temporarily REQUIRED]
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
    b = {};
    if ischar(varargin{1}) && strcmpi(varargin{1}, 'b')
        b = varargin(2:N+1);
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
    if numel(b) > 0 && numel(b) ~= N
        error('There should be as many b as intensity images')
    end
    if numel(bb) > 0 && numel(bb) ~= N
        error('There should be as many bounding boxes as intensity images')
    end

    % --- Parse inputs
    p = inputParser;
    p.FunctionName = 'updateMuLaplaceML';
    p.addParameter('lat',    [],    @isnumeric);
    p.addParameter('fwhm',   0,     @isnumeric);
    p.addParameter('loop',   '',    @(X) ischar(X) && any(strcmpi(X, {'slice', 'component', 'subject', 'none', ''})));
    p.addParameter('par',    false, @isscalar);
    p.addParameter('debug',  false, @isscalar);
    p.addParameter('output', []);
    p.parse(varargin{:});
    lat    = p.Results.lat;
    fwhm   = p.Results.fwhm;
    loop   = p.Results.loop;
    par    = p.Results.par;
    debug  = p.Results.debug;
    output = p.Results.output;
    
    if debug, fprintf('* updateMuLaplaceML\n'); end

    if isempty(lat)
        error('For now, the lattice size MUST be provided')
    end

    [par, loop] = autoParLoop(par, loop, isa(f{1}, 'file_array'), ...
                              size(f{1}, 3), size(f{1}, 4));
    if fwhm > 0 && strcmpi(loop, 'slice')
        loop = 'component';
    end
    
    switch lower(loop)
        case 'none'
            if debug, fprintf('   - No loop\n'); end
            mu = loopNone(f, c, b, bb, lat, output, fwhm);
        case 'component'
            if debug
                if par > 0
                    fprintf('   - Parallelise on components\n');
                else
                    fprintf('   - Serialise on components\n');
                end
            end
            mu = loopComponent(f, c, b, bb, lat, par, output, fwhm);
        case 'slice'
            if debug
                if par > 0
                    fprintf('   - Parallelise on slices\n');
                else
                    fprintf('   - Serialise on slices\n');
                end
            end
            mu = loopSlice(f, c, b, bb, lat, par, output);
        otherwise
            error('Unknown loop type ''%s''', loop)
    end

end

function mu = loopComponent(f, c, b, bb, lat, par, output, fwhm)

    nc = size(f{1}, 4);
    mu = prepareOnDisk(output, [lat nc], 'type', 'float32');
    
    if isempty(b)
        [b{1:numel(f)}] = deal(ones(1, nc));
    end
    if isempty(bb)
         bb1 = struct('x', 1:lat(1), 'y', 1:lat(2), 'z', 1:lat(3));
        [bb{1:numel(f)}] = deal(bb1);
    end

    if ~par
        for k=1:nc
            tmpf = zeros([lat numel(f)], 'single');
            tmpc = zeros([lat numel(f)], 'single');
            for n=1:numel(f)
                bx = bb{n}.x;
                by = bb{n}.y;
                bz = bb{n}.z;
                c1 = single(numeric(c{n}));
                tmpf(bx,by,bz,n) = single(f{n}(:,:,:,k)) ./ c1;
                tmpc(bx,by,bz,n) = c1./b{n}(k);
            end
            tmpf = smooth_gaussian(tmpf, fwhm);
            tmpc = smooth_gaussian(tmpc, fwhm);
            tmpf = wmedian(tmpf, tmpc);
            mu(:,:,:,k) = tmpf;
        end
    elseif isa(f{1}, 'file_array')
        parfor (k=1:nc, par)
            tmpf = zeros([lat numel(f)], 'single');
            tmpc = zeros([lat numel(f)], 'single');
            for n=1:numel(f)
                bx = bb{n}.x;
                by = bb{n}.y;
                bz = bb{n}.z;
                c1 = single(numeric(c{n}));
                tmpf(bx,by,bz,n) = single(slicevol(f{n}, k, 4)) ./ c1;
                tmpc(bx,by,bz,n) = c1./b{n}(k);
            end
            tmpf = smooth_gaussian(tmpf, fwhm);
            tmpc = smooth_gaussian(tmpc, fwhm);
            tmpf = wmedian(tmpf, tmpc);
            mu(:,:,:,k) = tmpf;
        end
    else
        parfor (k=1:nc, par)
            tmpf = zeros([lat numel(f)], 'single');
            tmpc = zeros([lat numel(f)], 'single');
            for n=1:numel(f)
                bx = bb{n}.x;
                by = bb{n}.y;
                bz = bb{n}.z;
                c1 = single(numeric(c{n}));
                tmpf(bx,by,bz,n) = single(f{n}(:,:,:,k)) ./ c1;
                tmpc(bx,by,bz,n) = c1./b{n}(k);
            end
            tmpf = smooth_gaussian(tmpf, fwhm);
            tmpc = smooth_gaussian(tmpc, fwhm);
            tmpf = wmedian(tmpf, tmpc);
            mu(:,:,:,k) = tmpf;
        end
    end
    
    if ~isempty(output)
        mu = saveOnDisk(output, mu, 'name', 'mu');
    end
end

function mu = loopSlice(f, c, b, bb, lat, par, output)

    nc = size(f{1}, 4);
    mu = prepareOnDisk(output, [lat nc], 'type', 'float32');
    
    if isempty(b)
        [b{1:numel(f)}] = deal(ones(1, nc));
    end
    if isempty(bb)
         bb1 = struct('x', 1:lat(1), 'y', 1:lat(2), 'z', 1:lat(3));
        [bb{1:numel(f)}] = deal(bb1);
    end
    
    % --- Compute mu
    if ~par
        for z=1:lat(3)
            tmpf = zeros([lat(1:2) 1 nc numel(f)], 'single');
            tmpc = zeros([lat(1:2) 1 nc numel(f)], 'single');
            for n=1:numel(f)
                bx = bb{n}.x;
                by = bb{n}.y;
                bz = bb{n}.z;
                fz = z-bz(1)+1;
                if fz >= 1 && fz <= size(f{n}, 3)
                    c1 = single(c{n}(:,:,fz));
                    tmpf(bx,by,1,:,n) = bsxfun(@rdivide, single(f{n}(:,:,fz,:)), c1);
                    tmpc(bx,by,1,:,n) = bsxfun(@rdivide, c1, reshape(b{n}, [1 1 1 nc]));
                end
            end
            tmpf = wmedian(tmpf, tmpc);
            mu(:,:,z,:) = tmpf;
        end
    elseif isa(f{1}, 'file_array')
        parfor (z=1:lat(3), par)
            tmpf = zeros([lat(1:2) 1 nc numel(f)], 'single');
            tmpc = zeros([lat(1:2) 1 nc numel(f)], 'single');
            for n=1:numel(f)
                bx = bb{n}.x;
                by = bb{n}.y;
                bz = bb{n}.z;
                fz = z-bz(1)+1;
                if fz >= 1 && fz <= size(f{n}, 3)
                    c1 = single(slicevol(c{n}, fz, 3));
                    tmpf(bx,by,1,:,n) = bsxfun(@rdivide, single(slicevol(f{n}, fz, 3)), c1);
                    tmpc(bx,by,1,:,n) = bsxfun(@rdivide, c1, reshape(b{n}, [1 1 1 nc]));
                end
            end
            tmpf = wmedian(tmpf, tmpc);
            mu(:,:,z,:) = tmpf;
        end
    else
        parfor (z=1:lat(3), par)
            tmpf = zeros([lat(1:2) 1 nc numel(f)], 'single');
            tmpc = zeros([lat(1:2) 1 nc numel(f)], 'single');
            for n=1:numel(f)
                bx = bb{n}.x;
                by = bb{n}.y;
                bz = bb{n}.z;
                fz = z-bz(1)+1;
                if fz >= 1 && fz <= size(f{n}, 3)
                    c1 = single(c{n}(:,:,fz));
                    tmpf(bx,by,1,:,n) = bsxfun(@rdivide, single(f{n}(:,:,fz,:)), c1);
                    tmpc(bx,by,1,:,n) = bsxfun(@rdivide, c1, reshape(b{n}, [1 1 1 nc]));
                    tmpf(:,:,1,:,n) = single(f{n}(:,:,fz,:));
                end
            end
            tmpf = wmedian(tmpf, tmpc);
            mu(:,:,z,:) = tmpf;
        end
    end
    
    if ~isempty(output)
        mu = saveOnDisk(output, mu, 'name', 'mu');
    end
end

function mu = loopNone(f, c, b, bb, lat, output, fwhm)

    if nargin < 7
        fwhm = 0;
    end

    nc = size(f{1}, 4);
    
    if isempty(b)
        [b{1:numel(f)}] = deal(ones(1, nc));
    end
    if isempty(bb)
         bb1 = struct('x', 1:lat(1), 'y', 1:lat(2), 'z', 1:lat(3));
        [bb{1:numel(f)}] = deal(bb1);
    end

    mu   = zeros([lat nc numel(f)], 'single');
    tmpc = zeros([lat nc numel(f)], 'single');
    for n=1:numel(f)
        bx = bb{n}.x;
        by = bb{n}.y;
        bz = bb{n}.z;
        c1 = single(numeric(c{n}));
        mu(bx,by,bz,:,n)   = bsxfun(@rdivide, single(numeric(f{n})), c1);
        tmpc(bx,by,bz,:,n) = bsxfun(@rdivide, c1, reshape(b{n}, [1 1 1 nc]));
    end
    mu   = smooth_gaussian(mu, fwhm);
    tmpc = smooth_gaussian(tmpc, fwhm);
    mu   = wmedian(mu, tmpc);
    
    if ~isempty(output)
        mu = saveOnDisk(output, mu, 'name', 'mu');
    end
end

function m = wmedian(f, c)
% Weighted median along last non-singleton dimension

    sizef = size(f);
    dimc  = length(size(c));

    [f, ind] = sort(f,dimc);
    c = c(ind);
    c = bsxfun(@rdivide, cumsum(c, dimc), sum(c, dimc)) > 0.5;
    [~, ind] = max(c, [], dimc);
        
    f = reshape(f, [], sizef(dimc));
    f = f(sub2ind(size(f), 1:size(f,1), ind(:)'));
    f = reshape(f, sizef(1:end-1));

    m = f;
end