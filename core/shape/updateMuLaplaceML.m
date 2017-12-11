function mu = updateMuLaplaceML(varargin)
% FORMAT mu = updateMuLaplaceML('f', pf1, ..., pfN, 
%                               'c', c1, ..., cN,
%                               ('b', b1, ..., bN),
%                               ('fwhm', fwhm),
%                               ('loop', loop), ('par', par))
% ** Required **
% pfi  - Image pushed in template space
% ci   - Pushed voxel count
% ** Optional **
% bi   - Individual noise variance (b)
% ** Keyword arguments **
% fwhm  - Smoothing kernel used as pseudo prior [do not use]
% loop - How to split processing: 'component', 'slice', 'none' or '' [auto]
% par  - Distribute compute [auto]
% ** Output **
% mu   - Updated template
%
% Closed form M-step update of the template for the Laplace noise model.
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
    if ischar(varargin{1}) && strcmpi(varargin{1}, 'b')
        b = varargin(2:N+1);
        varargin = varargin(N+2:end);
    end

    if numel(c) ~= N
        error('There should be as many count as intensity images')
    end
    if numel(b) > 0 && numel(b) ~= N
        error('There should be as many b as intensity images')
    end

    % --- Parse inputs
    p = inputParser;
    p.FunctionName = 'updateMuLaplaceML';
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
    
    if debug, fprintf('* updateMuLaplaceML\n'); end;

    [par, loop] = autoParLoop(par, loop, isa(f{1}, 'file_array'), ...
                              size(f{1}, 3), size(f{1}, 4));
    if fwhm > 0 && strcmpi(loop, 'slice')
        loop = 'component';
    end
    
    switch lower(loop)
        case 'none'
            if debug, fprintf('   - No loop\n'); end;
            mu = loopNone(f, c, output, fwhm);
        case 'component'
            if debug
                if par > 0
                    fprintf('   - Parallelise on components\n');
                else
                    fprintf('   - Serialise on components\n');
                end
            end
            mu = loopComponent(f, c, par, output, fwhm);
        case 'slice'
            if debug
                if par > 0
                    fprintf('   - Parallelise on slices\n');
                else
                    fprintf('   - Serialise on slices\n');
                end
            end
            mu = loopSlice(f, c, par, output);
        otherwise
            error('Unknown loop type ''%s''', loop)
    end

end

function mu = loopComponent(f, c, par, output, fwhm)

    mu = prepareOnDisk(output, size(f{1}), 'type', 'float32');
    dim = [size(mu) 1 1];
    lat = dim(1:3);
    nc  = dim(4);
    
    if ~par
        for k=1:nc
            tmpf = zeros([lat numel(f)], 'single');
            tmpc = zeros([lat numel(f)], 'single');
            for n=1:numel(f)
                tmpf(:,:,:,n) = single(f{n}(:,:,:,k));
                tmpc(:,:,:,n) = single(numeric(c{n}));
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
                tmpf(:,:,:,n) = single(slicevol(f{n}, k, 4));
                tmpc(:,:,:,n) = single(numeric(c{n}));
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
                tmpf(:,:,:,n) = single(f{n}(:,:,:,k));
                tmpc(:,:,:,n) = single(numeric(c{n}));
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

function mu = loopSlice(f, c, par, output)

    mu = prepareOnDisk(output, size(f{1}), 'type', 'float32');
    dim = [size(mu) 1 1];
    lat = dim(1:3);
    nc  = dim(4);
    
    % --- Compute count
    tmpc = zeros(lat, 'single');
    for n=1:numel(f)
        tmpc = tmpc + single(numeric(c{n}));
    end
    
    % --- Compute mu
    if ~par
        for z=1:dim(3)
            tmpf = zeros([lat(1:2) 1 nc numel(f)], 'single');
            for n=1:numel(f)
                tmpf(:,:,1,:,n) = single(f{n}(:,:,z,:));
            end
            tmpf = wmedian(tmpf, tmpc(:,:,z));
            mu(:,:,z,:) = tmpf;
        end
    elseif isa(f{1}, 'file_array')
        parfor (z=1:dim(3), par)
            tmpf = zeros([lat(1:2) 1 nc numel(f)], 'single');
            for n=1:numel(f)
                tmpf(:,:,1,:,n) = single(slicevol(f{n}, z, 3));
            end
            tmpf = wmedian(tmpf, tmpc(:,:,z));
            mu(:,:,z,:) = tmpf;
        end
    else
        parfor (z=1:dim(3), par)
            tmpf = zeros([lat(1:2) 1 nc numel(f)], 'single');
            for n=1:numel(f)
                tmpf(:,:,1,:,n) = single(f{n}(:,:,z,:));
            end
            tmpf = wmedian(tmpf, tmpc(:,:,z));
            mu(:,:,z,:) = tmpf;
        end
    end
    
    if ~isempty(output)
        mu = saveOnDisk(output, mu, 'name', 'mu');
    end
end

function mu = loopNone(f, c, output)

    dim = [size(f{1}) 1 1];
    lat = dim(1:3);
    nc  = dim(4);
    
    mu   = zeros([lat nc numel(f)], 'single');
    tmpc = zeros([lat numel(f)], 'single');
    for n=1:numel(f)
        mu(:,:,:,:,n) = single(numeric(f{n}));
        tmpc(:,:,:,n) = single(numeric(c{n}));
    end
    mu   = smooth_gaussian(mu, fwhm);
    tpmc = smooth_gaussian(tmpc, fwhm);
    mu   = wmedian(mu, tmpc);
    
    if ~isempty(output)
        mu = saveOnDisk(output, mu, 'name', 'mu');
    end
end

function m = wmedian(f, c)
% Weighted median along last non-singleton dimension

    sizef = size(f);
    dimc  = length(size(c));

    [f, ind] = sort(f, sizef(end));
    c = c(:,:,:,ind);
    c = bsxfun(@rdivide, cumsum(c, dimc), sum(c, dimc)) > 0.5;
    [~, ind] = max(c, [], 4);
        
    f = reshape(f, [], sizef(end));
    f = f(sub2ind(size(f), 1:size(f,1), ind));
    f = reshape(f, sizef(1:end-1));

    m = f;
end