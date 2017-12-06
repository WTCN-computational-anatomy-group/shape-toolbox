function mu = updateMuBernoulliML(varargin)
% FORMAT mu = updateMuBernoulliML(pf1, ..., pfN, c1, ..., cN,
%                                 ('fwhm', fwhm),
%                                 ('loop', loop), ('par', par))
% ** Required **
% pfi  - Image pushed in template space
% ci   - Pushed voxel count
% ** Keyword arguments **
% fwhm  - Smoothing kernel used as pseudo prior [do not use]
% loop - How to split processing: 'slice', 'none' or '' [auto]
% par  - Distribute compute [auto]
% ** Output **
% mu   - Updated template
%
% Closed form M-step update of the template for the Bernoulli model.
% (Maximum likelihood update: no prior on Mu)

    N = 0;
    while N < numel(varargin) && ~ischar(varargin{N+1})
        N = N+1;
    end
    if mod(N,2)
        error('There should be as many intensity as count images')
    end
    f = varargin(1:(N/2));
    c = varargin((N/2+1):N);
    varargin = varargin(N+1:end);

    % --- Parse inputs
    p = inputParser;
    p.FunctionName = 'updateMuBernoulliML';
    p.addParameter('fwhm',   0,     @isnumeric);
    p.addParameter('loop',   '',    @(X) ischar(X) && any(strcmpi(X, {'slice', 'subject', 'none', ''})));
    p.addParameter('par',    false, @isscalar);
    p.addParameter('debug',  false, @isscalar);
    p.addParameter('output', false);
    p.parse(varargin{:});
    fwhm   = p.Results.fwhm;
    loop   = p.Results.loop;
    par    = p.Results.par;
    debug  = p.Results.debug;
    
    if debug, fprintf('* updateMuBernoulliML\n'); end;

    [par, loop] = autoParLoop(par, loop, isa(f{1}, 'file_array'), ...
                              size(f{1}, 3), size(f{1}, 4));
    if fwhm > 0 && strcmpi(loop, 'slice')
        loop = 'none';
    end
    
    switch lower(loop)
        case 'none'
            if debug, fprintf('   - No loop\n'); end;
            mu = loopNone(f, c, output, fwhm);
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

function mu = loopSlice(f, c, par, output)

    mu = prepareOnDisk(output, size(f{1}), 'type', 'float32');
    dim = [size(mu) 1];
    lat = dim(1:3);
    
    % --- Compute count
    tmpc = zeros(lat, 'single');
    for n=1:numel(f)
        tmpc = tmpc + single(numeric(c{n}));
    end
    
    % -- Compute mu
    parfor (z=1:dim(3), par)
        tmpf = zeros(lat(1:2), 'single');
        for n=1:numel(f)
            tmpf = tmpf + single(f{n}(:,:,z,:));
        end
        tmpf = tmpf ./ tmpc(:,:,z);
        tmpf = max(1-eps('single'), min(eps('single'), tmpf));
        tmpf(~isfinite(tmpf)) = 0.5;
        tmpf = log(tmpf) - log(1 - tmpf);
        mu(:,:,z,:) = tmpf;
    end
    
    if ~isempty(output)
        mu = saveOnDisk(output, mu, 'name', 'mu');
    end
end

function mu = loopNone(f, c, output)

    dim = [size(f{1}) 1];
    lat = dim(1:3);
    
    mu   = zeros(lat, 'single');
    tmpc = zeros(lat, 'single');
    for n=1:numel(f)
        mu   = mu   + single(numeric(f{n}));
        tmpc = tmpc + single(numerci(c{n}));
    end
    mu   = smooth_gaussian(mu, fwhm);
    tmpc = smooth_gaussian(tmpc, fwhm);
    mu = mu ./ tmpc;
    mu = max(1-eps('single'), min(eps('single'), mu));
    mu(~isfinite(mu)) = 0.5;
    mu = log(mu) - log(1 - mu);
    
    if ~isempty(output)
        mu = saveOnDisk(output, mu, 'name', 'mu');
    end
end