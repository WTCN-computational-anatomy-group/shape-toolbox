function mu = updateMuBernoulliML(varargin)
%__________________________________________________________________________
%
% Closed form M-step update of the template for the Bernoulli noise model.
% (Maximum likelihood update: no prior on Mu)
%
% -------------------------------------------------------------------------
%
% FORMAT mu = updateMuNormalML('f', ..., 'c', ..., ('bb', ...), ...)
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
    p.FunctionName = 'updateMuBernoulliML';
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
    output = p.Results.output;
    debug  = p.Results.debug;
    
    if debug, fprintf('* updateMuBernoulliML\n'); end

    if isempty(lat)
        error('For now, the lattice size MUST be provided')
    end

    [par, loop] = autoParLoop(par, loop, isa(f{1}, 'file_array'), ...
                              size(f{1}, 3), size(f{1}, 4));
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

    mu = prepareOnDisk(output, lat, 'type', 'float32');
    
    if isempty(bb)
         bb1 = struct('x', 1:lat(1), 'y', 1:lat(2), 'z', 1:lat(3));
        [bb{1:numel(f)}] = deal(bb1);
    end

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
    
    % -- Compute mu
    if ~par
        for z=1:lat(3)
            tmpf = zeros(lat(1:2), 'single');
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
                    tmpf(bx,by) = tmpf(bx,by) + single(f{n}(:,:,fz));
                end
            end
            tmpf = tmpf ./ tmpc(:,:,z);
            tmpf = min(1-eps('single'), max(eps('single'), tmpf));
            tmpf(~isfinite(tmpf)) = 0.5;
            tmpf = log(tmpf) - log(1 - tmpf);
            mu(:,:,z,:) = tmpf;
        end
    elseif isa(f{1}, 'file_array')
        parfor (z=1:lat(3), par)
            tmpf = zeros(lat(1:2), 'single');
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
                    tmpf(bx,by) = tmpf(bx,by) + single(slice(f{n}, fz, 3));
                end
            end
            tmpf = tmpf ./ tmpc(:,:,z);
            tmpf = min(1-eps('single'), max(eps('single'), tmpf));
            tmpf(~isfinite(tmpf)) = 0.5;
            tmpf = log(tmpf) - log(1 - tmpf);
            mu(:,:,z,:) = tmpf;
        end
    else
        parfor (z=1:lat(3), par)
            tmpf = zeros(lat(1:2), 'single');
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
                    tmpf(bx,by) = tmpf(bx,by) + single(f{n}(:,:,fz));
                end
            end
            tmpf = tmpf ./ tmpc(:,:,z);
            tmpf = min(1-eps('single'), max(eps('single'), tmpf));
            tmpf(~isfinite(tmpf)) = 0.5;
            tmpf = log(tmpf) - log(1 - tmpf);
            mu(:,:,z,:) = tmpf;
        end
    end
    
    if ~isempty(output)
        mu = saveOnDisk(output, mu, 'name', 'mu');
    end
end

function mu = loopNone(f, c, bb, output, fwhm)
    
    if isempty(bb)
         bb1 = struct('x', 1:lat(1), 'y', 1:lat(2), 'z', 1:lat(3));
        [bb{1:numel(f)}] = deal(bb1);
    end

    mu   = zeros(lat, 'single');
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
        mu(bx,by,bz)   = mu(bx,by,bz)   + single(numeric(f{n}));
        tmpc(bx,by,bz) = tmpc(bx,by,bz) + single(numeric(c{n}));
    end
    mu   = smooth_gaussian(mu, fwhm);
    tmpc = smooth_gaussian(tmpc, fwhm);
    mu = mu ./ tmpc;
    mu = min(1-eps('single'), mimaxn(eps('single'), mu));
    mu(~isfinite(mu)) = 0.5;
    mu = log(mu) - log(1 - mu);
    
    if ~isempty(output)
        mu = saveOnDisk(output, mu, 'name', 'mu');
    end
end