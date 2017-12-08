function mu = updateMuNormalML2(varargin)
% FORMAT mu = updateMuNormalML('f', pf1, ..., pfN, 
%                              'c', c1, ..., cN,
%                              ('s', s1, ..., sN), 
%                              ('loop', loop), ('par', par))
% ** Required **
% pfi  - Image pushed in template space
% ci   - Pushed voxel count
% ** Optional **
% si   - Individual noise variance (sigma2)
% ** Keyword arguments **
% fwhm - Smoothing kernel used as pseudo prior [do not use]
% loop - How to split processing: 'component', 'slice', 'none' or '' [auto]
% par  - Distribute compute [auto]
% ** Output **
% mu   - Updated template
%
% Closed form M-step update of the template for the Gaussian noise model.
% (Maximum likelihood update: no prior on Mu)


    i = 1;
    if ischar(varargin{1})
        i = i + 1;
        N = 0;
        while i <= numel(varargin) && ~ischar(varargin{i})
            i = i + 1;
            N = N + 1;
        end
    end
    
    if ischar(varargin{1}) && strcmpi(varargin{1}, 'f')
        f = varargin(2:N+1);
        varargin = varargin(N+2:end);
    end
    if ischar(varargin{1}) && strcmpi(varargin{1}, 'c')
        c = varargin(2:N+1);
        varargin = varargin(N+2:end);
    end
    if ischar(varargin{1}) && strcmpi(varargin{1}, 's')
        s = varargin(2:N+1);
        varargin = varargin(N+2:end);
    end

    if numel(c) ~= N
        error('There should be as many count as intensity images')
    end
    if numel(s) > 0 && numel(s) ~= N
        error('There should be as many sigma2 as intensity images')
    end

    % --- Parse inputs
    p = inputParser;
    p.FunctionName = 'updateMuNormalML';
    p.addParameter('fwhm',   0,     @isnumeric);
    p.addParameter('loop',   '',    @(X) ischar(X) && any(strcmpi(X, {'slice', 'component', 'subject', 'none', ''})));
    p.addParameter('par',    false, @isscalar);
    p.addParameter('debug',  false, @isscalar);
    p.addParameter('output', false);
    p.parse(varargin{:});
    fwhm   = p.Results.fwhm;
    loop   = p.Results.loop;
    par    = p.Results.par;
    debug  = p.Results.debug;
    output = p.Results.output;
    
    if debug, fprintf('* updateMuNormalML\n'); end;

    [par, loop] = autoParLoop(par, loop, isa(f{1}, 'file_array'), ...
                              size(f{1}, 3), size(f{1}, 4));
    if fwhm > 0 && strcmpi(loop, 'slice')
        loop = 'component';
    end
    
    switch lower(loop)
        case 'none'
            if debug, fprintf('   - No loop\n'); end;
            mu = loopNone(f, c, s, output, fwhm);
        case 'component'
            if debug
                if par > 0
                    fprintf('   - Parallelise on components\n');
                else
                    fprintf('   - Serialise on components\n');
                end
            end
            mu = loopComponent(f, c, s, par, output, fwhm);
        case 'slice'
            if debug
                if par > 0
                    fprintf('   - Parallelise on slices\n');
                else
                    fprintf('   - Serialise on slices\n');
                end
            end
            mu = loopSlice(f, c, s, par, output);
        otherwise
            error('Unknown loop type ''%s''', loop)
    end

end

function mu = loopComponent(f, c, s, par, output, fwhm)

    mu = prepareOnDisk(output, size(f{1}), 'type', 'float32');
    dim = [size(mu) 1 1];
    lat = dim(1:3);
    nc  = dim(4);
    
    tmp = struct('f', f, 'c', c);
    if ~isempty(s)
        [tmp.s] = deal(s{:});
    end
    if ~par
        for k=1:nc, par
            tmpf = zeros(lat, 'single');
            tmpc = zeros(lat, 'single');
            for n=1:numel(tmp)
                if ~isempty(s)
                    s1 = tmp(n).s(k);
                    c1 = single(numeric(tmp(n).c)/s1);
                else
                    c1 = single(numeric(tmp(n).c));
                end
                tmpc = tmpc + c1;
                tmpf = tmpf + c1 .* single(tmp(n).f(:,:,:,k));
            end
            tmpf = smooth_gaussian(tmpf, fwhm);
            tmpc = smooth_gaussian(tmpc, fwhm);
            tmpf = tmpf ./ tmpc;
            mu(:,:,:,k) = tmpf;
        end
    elseif isa(tmp(1).f, 'file_array')
        parfor (k=1:nc, par)
            tmpf = zeros(lat, 'single');
            tmpc = zeros(lat, 'single');
            for n=1:numel(tmp)
                if ~isempty(s)
                    s1 = tmp(n).s(k);
                    c1 = single(numeric(tmp(n).c)/s1);
                else
                    c1 = single(numeric(tmp(n).c));
                end
                tmpc = tmpc + c1;
                tmpf = tmpf + c1 .* single(slicevol(tmp(n).f, k, 4));
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
                if ~isempty(s)
                    s1 = tmp(n).s(k);
                    c1 = single(numeric(tmp(n).c)/s1);
                else
                    c1 = single(numeric(tmp(n).c));
                end
                tmpc = tmpc + c1;
                tmpf = tmpf + c1 .* single(tmp(n).f(:,:,:,k));
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

function mu = loopSlice(f, c, s, par, output)

    mu = prepareOnDisk(output, size(f{1}), 'type', 'float32');
    dim = [size(mu) 1 1];
    lat = dim(1:3);
    nc  = dim(4);
    
    % --- Compute mu
    tmp = struct('f', f, 'c', c);
    if ~isempty(s)
        [tmp.s] = deal(s{:});
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
                if isempty(s)
                    c1 = single(numeric(tmp(n).c(:,:,z)));
                    f1 = single(tmp(n).f(:,:,z,:));
                else
                    s1 = reshape(tmp(n).s, [1 1 1 size(f1, 4)]);
                    c1 = bsxfun(@rdivide, ...
                                single(numeric(tmp(n).c(:,:,z))), s1);
                    f1 = bsxfun(@times, single(tmp(n).f(:,:,z,:)), s1);
                end
                tmpc = tmpc + c1;
                tmpf = tmpf + f1;
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
                if isempty(s)
                    c1 = single(numeric(slicevol(tmp(n).c, z, 3)));
                    f1 = single(slicevol(tmp(n).f, z, 3));
                else
                    s1 = reshape(tmp(n).s, [1 1 1 size(f1, 4)]);
                    c1 = bsxfun(@rdivide, ...
                                single(numeric(slicevol(tmp(n).c, z, 3))), s1);
                    f1 = bsxfun(@times, single(slicevol(tmp(n).f, z, 3)), s1);
                end
                tmpc = tmpc + c1;
                tmpf = tmpf + f1;
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
                if isempty(s)
                    c1 = single(numeric(tmp(n).c(:,:,z)));
                    f1 = single(tmp(n).f(:,:,z,:));
                else
                    s1 = reshape(tmp(n).s, [1 1 1 size(f1, 4)]);
                    c1 = bsxfun(@rdivide, ...
                                single(numeric(tmp(n).c(:,:,z))), s1);
                    f1 = bsxfun(@times, single(tmp(n).f(:,:,z,:)), s1);
                end
                tmpc = tmpc + c1;
                tmpf = tmpf + f1;
            end
            mu(:,:,z,:) = bsxfun(@rdivide, tmpf, tmpc);
        end
    end
    
    if ~isempty(output)
        mu = saveOnDisk(output, mu, 'name', 'mu');
    end
end

function mu = loopNone(f, c, s, output, fwhm)

    dim = [size(f{1}) 1 1];
    lat = dim(1:3);
    nc  = dim(4);
    
    mu   = zeros([lat nc], 'single');
    if isempty(s)
        tmpc = zeros(lat, 'single');
    else
        tmpc = zeros([lat nc], 'single');
    end
    for n=1:numel(f)
        if ~isempty(s)
            s1 = reshape(s{n}, [1 1 1 nc]);
            c1 = bsxfun(@rdivide, single(numeric(c{n})), s1);
            f1   = bsxfun(@times, single(numeric(f{n})), s1);
        else
            c1 = single(numeric(c{n}));
            f1 = single(numeric(f{n}));
        end
        tmpc = tmpc + c1;
        mu   = mu + f1;
    end
    mu   = smooth_gaussian(mu, fwhm);
    tmpc = smooth_gaussian(tmpc, fwhm);
    mu   = bsxfun(@rdivide, mu, tmpc);
    
    if ~isempty(output)
        mu = saveOnDisk(output, mu, 'name', 'mu');
    end
end