function mu = updateMuNormalML(varargin)
%__________________________________________________________________________
%
% Closed form M-step update of the template for the Gaussian noise model.
% (Maximum likelihood update: no prior on Mu)
%
% -------------------------------------------------------------------------
%
% FORMAT mu = updateMuNormalML('f', ..., 'c', ...,
%                              ('s', ...),  ('bb', ...), ...)
%
% REQUIRED KEYWORD (LIST)
% -----------------------
% f   - Images pushed in template space
% c   - Pushed voxel counts
%
% OPTIONAL KEYWORD (LIST)
% -----------------------
% s    - Individual noise variances (sigma2)
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
    s = {};
    if ischar(varargin{1}) && strcmpi(varargin{1}, 's')
        s = varargin(2:N+1);
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
    if numel(s) > 0 && numel(s) ~= N
        error('There should be as many sigma2 as intensity images')
    end
    if numel(bb) > 0 && numel(bb) ~= N
        error('There should be as many bounding boxes as intensity images')
    end


    % --- Parse inputs
    p = inputParser;
    p.FunctionName = 'updateMuNormalML';
    p.addParameter('lat',    [],    @isnumeric);
    p.addParameter('fwhm',   0,     @isnumeric);
    p.addParameter('loop',   '',    @(X) ischar(X) && any(strcmpi(X, {'slice', 'component', 'subject', 'none', ''})));
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
    
    if debug, fprintf('* updateMuNormalML\n'); end;

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
            if debug, fprintf('   - No loop\n'); end;
            mu = loopNone(f, c, s, bb, lat, output, fwhm);
        case 'component'
            if debug
                if par > 0
                    fprintf('   - Parallelise on components\n');
                else
                    fprintf('   - Serialise on components\n');
                end
            end
            mu = loopComponent(f, c, s, bb, lat, par, output, fwhm);
        case 'slice'
            if debug
                if par > 0
                    fprintf('   - Parallelise on slices\n');
                else
                    fprintf('   - Serialise on slices\n');
                end
            end
            mu = loopSlice(f, c, s, bb, lat, par, output);
        otherwise
            error('Unknown loop type ''%s''', loop)
    end

end

function mu = loopComponent(f, c, s, bb, lat, par, output, fwhm)

    nc  = size(f{1}, 4);
    mu  = prepareOnDisk(output, [lat nc], 'type', 'float32');
    
    tmp = struct('f', f, 'c', c);
    if ~isempty(s)
        [tmp.s] = deal(s{:});
    else
        [tmp.s] = deal(ones(nc, 1));
    end
    if isempty(bb)
        bb = struct('x', 1:lat(1), 'y', 1:lat(2), 'z', 1:lat(3));
        [tmp.bb] = deal(bb);
    else
        [tmp.bb] = deal(bb{:});
    end
    if ~par
        for k=1:nc
            tmpf = zeros(lat, 'single');
            tmpc = zeros(lat, 'single');
            for n=1:numel(tmp)
                s1 = tmp(n).s(k);
                c1 = single(numeric(tmp(n).c)/s1);
                bx = tmp(n).bb.x;
                by = tmp(n).bb.y;
                bz = tmp(n).bb.z;
                tmpc(bx,by,bz) = tmpc(bx,by,bz) + c1;
                tmpf(bx,by,bz) = tmpf(bx,by,bz) + c1 .* single(f{n}(:,:,:,k));
            end
            tmpf = smooth_gaussian(tmpf, fwhm);
            tmpc = smooth_gaussian(tmpc, fwhm);
            tmpf = tmpf ./ tmpc;
            mu(tmp(n).bb.x,tmp(n).bb.y,tmp(n).bb.z,k) = tmpf;
        end
    elseif isa(tmp(1).f, 'file_array')
        parfor (k=1:nc, par)
            tmpf = zeros(lat, 'single');
            tmpc = zeros(lat, 'single');
            for n=1:numel(tmp)
                s1 = tmp(n).s(k);
                c1 = single(numeric(tmp(n).c)/s1);
                bx = tmp(n).bb.x;
                by = tmp(n).bb.y;
                bz = tmp(n).bb.z;
                tmpc(bx,by,bz) = tmpc(bx,by,bz) + c1;
                tmpf(bx,by,bz) = tmpf(bx,by,bz) + c1 .* single(slicevol(tmp(n).f, k, 4));
            end
            tmpf = smooth_gaussian(tmpf, fwhm);
            tmpc = smooth_gaussian(tmpc, fwhm);
            tmpf = tmpf ./ tmpc;
            mu(:,:,:,k) = tmpf;
        end
    else
        parfor (k=1:nc, par)
            tmpf = zeros(lat, 'single');
            tmpc = zeros(lat, 'single');
            for n=1:numel(tmp)
                s1 = tmp(n).s(k);
                c1 = single(numeric(tmp(n).c)/s1);
                bx = tmp(n).bb.x;
                by = tmp(n).bb.y;
                bz = tmp(n).bb.z;
                tmpc(bx,by,bz) = tmpc(bx,by,bz) + c1;
                tmpf(bx,by,bz) = tmpf(bx,by,bz) + c1 .* single(tmp(n).f(:,:,:,k));
            end
            tmpf = smooth_gaussian(tmpf, fwhm);
            tmpc = smooth_gaussian(tmpc, fwhm);
            tmpf = tmpf ./ tmpc;
            mu(:,:,:,k) = tmpf;
        end
    end
    
    if ~isempty(output)
        mu = saveOnDisk(output, mu, 'name', 'mu');
    end
end

function mu = loopSlice(f, c, s, bb, lat, par, output)

    nc  = size(f{1}, 4);
    mu = prepareOnDisk(output, [lat nc], 'type', 'float32');
    
    % --- Compute mu
    tmp = struct('f', f, 'c', c);
    if ~isempty(s)
        [tmp.s] = deal(s{:});
    else
        [tmp.s] = deal(ones(nc, 1));
    end
    if isempty(bb)
        bb = struct('x', 1:lat(1), 'y', 1:lat(2), 'z', 1:lat(3));
        [tmp.bb] = deal(bb);
    else
        [tmp.bb] = deal(bb{:});
    end
    if ~par
        for z=1:dim(3)
            tmpf = zeros([lat(1:2) 1 nc], 'single');
            if isempty(s)
                tmpc = zeros([lat(1:2) 1],'single');
            else
                tmpc = zeros([lat(1:2) 1 nc],'single');
            end
            for n=1:numel(tmp)
                bx = tmp(n).bb.x;
                by = tmp(n).bb.y;
                bz = tmp(n).bb.z;
                fz = z-bz(1)+1;
                if fz >= 1 && fz <= size(tmp(1).f, 3)
                    s1 = reshape(tmp(n).s, [1 1 1 size(tmp(n).f, 4)]);
                    c1 = bsxfun(@rdivide, ...
                                single(numeric(tmp(n).c(:,:,fz))), s1);
                    f1 = bsxfun(@times, single(tmp(n).f(:,:,fz,:)), s1);
                    tmpc(bx,by,z) = tmpc(bx,by,z) + c1;
                    tmpf(bx,by,z) = tmpf(bx,by,z) + f1;
                end
            end
            mu(:,:,z,:) = bsxfun(@rdivide, tmpf, tmpc);
        end
    elseif isa(tmp(1).f, 'file_array')
        parfor (z=1:dim(3), par)
            tmpf = zeros([lat(1:2) 1 nc], 'single');
            if isempty(s)
                tmpc = zeros([lat(1:2) 1],'single');
            else
                tmpc = zeros([lat(1:2) 1 nc],'single');
            end
            for n=1:numel(tmp)
                bx = tmp(n).bb.x;
                by = tmp(n).bb.y;
                bz = tmp(n).bb.z;
                fz = z-bz(1)+1;
                if fz >= 1 && fz <= size(tmp(n).f, 3)
                    s1 = reshape(tmp(n).s, [1 1 1 size(tmp(n).f, 4)]);
                    c1 = bsxfun(@rdivide, ...
                                single(numeric(slicevol(tmp(n).c, fz, 3))), s1);
                    f1 = bsxfun(@times, single(slicevol(tmp(n).f, fz, 3)), s1);
                    tmpc(bx,by,z) = tmpc(bx,by,z) + c1;
                    tmpf(bx,by,z) = tmpf(bx,by,z) + f1;
                end
            end
            mu(:,:,z,:) = bsxfun(@rdivide, tmpf, tmpc);
        end
    else
        parfor (z=1:dim(3), par)
            tmpf = zeros([lat(1:2) 1 nc], 'single');
            if isempty(s)
                tmpc = zeros([lat(1:2) 1],'single');
            else
                tmpc = zeros([lat(1:2) 1 nc],'single');
            end
            for n=1:numel(tmp)
                bx = tmp(n).bb.x;
                by = tmp(n).bb.y;
                bz = tmp(n).bb.z;
                fz = z-bz(1)+1;
                if fz >= 1 && fz <= size(tmp(n).f, 3)
                    s1 = reshape(tmp(n).s, [1 1 1 size(tmp(n).f, 4)]);
                    c1 = bsxfun(@rdivide, ...
                                single(numeric(tmp(n).c(:,:,fz))), s1);
                    f1 = bsxfun(@times, single(tmp(n).f(:,:,fz,:)), s1);
                    tmpc(bx,by,z) = tmpc(bx,by,z) + c1;
                    tmpf(bx,by,z) = tmpf(bx,by,z) + f1;
                end
            end
            mu(:,:,z,:) = bsxfun(@rdivide, tmpf, tmpc);
        end
    end
    
    if ~isempty(output)
        mu = saveOnDisk(output, mu, 'name', 'mu');
    end
end

function mu = loopNone(f, c, s, bb, lat, output, fwhm)

    nc  = size(f{1}, 4);
    
    mu   = zeros([lat nc], 'single');
    if isempty(s)
        tmpc = zeros(lat, 'single');
    else
        tmpc = zeros([lat nc], 'single');
    end
    for n=1:numel(f)
        bx = bb{n}.x;
        by = bb{n}.y;
        bz = bb{n}.z;
        if ~isempty(s)
            s1 = reshape(s{n}, [1 1 1 nc]);
            c1 = bsxfun(@rdivide, single(numeric(c{n})), s1);
            f1 = bsxfun(@times, single(numeric(f{n})), s1);
        else
            c1 = single(numeric(c{n}));
            f1 = single(numeric(f{n}));
        end
        tmpc(bx,by,bz,:) = tmpc(bx,by,bz,:) + c1;
        mu(bx,by,bz,:)   = mu(bx,by,bz,:)   + f1;
    end
    mu   = smooth_gaussian(mu, fwhm);
    tmpc = smooth_gaussian(tmpc, fwhm);
    mu   = bsxfun(@rdivide, mu, tmpc);
    
    if ~isempty(output)
        mu = saveOnDisk(output, mu, 'name', 'mu');
    end
end