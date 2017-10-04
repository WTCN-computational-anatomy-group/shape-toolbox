function mu = updateMuNormalML(varargin)
% FORMAT mu = updateMuNormalML(pf1, ..., pfN, c1, ..., cN,
%                              ('loop', loop), ('par', par))
% ** Required **
% pfi  - Image pushed in template space
% ci   - Pushed voxel count
% ** Keyword arguments **
% loop - How to split processing: 'component', 'slice', 'none' or '' [auto]
% par  - Distribute compute [auto]
% ** Output **
% mu   - Updated template
%
% Closed form M-step update of the template for the Gaussian noise model.
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
    p.FunctionName = 'updateMuNormalML';
    p.addParameter('loop',   '',    @(X) ischar(X) && any(strcmpi(X, {'slice', 'component', 'subject', 'none', ''})));
    p.addParameter('par',    false, @isscalar);
    p.addParameter('debug',  false, @isscalar);
    p.addParameter('output', false);
    p.parse(varargin{:});
    loop   = p.Results.loop;
    par    = p.Results.par;
    debug  = p.Results.debug;
    output = p.Results.output;
    
    if debug, fprintf('* updateMuNormalML\n'); end;

    [par, loop] = autoParLoop(par, loop, isa(f{1}, 'file_array'), ...
                              size(f{1}, 3), size(f{1}, 4));
    
    switch lower(loop)
        case 'none'
            if debug, fprintf('   - No loop\n'); end;
            mu = loopNone(f, c, output);
        case 'component'
            if debug
                if par > 0
                    fprintf('   - Parallelise on components\n');
                else
                    fprintf('   - Serialise on components\n');
                end
            end
            mu = loopComponent(f, c, par, output);
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

function mu = loopComponent(f, c, par, output)

    mu = prepareOnDisk(output, size(f{1}), 'type', 'float32');
    dim = [size(mu) 1 1];
    lat = dim(1:3);
    nc  = dim(4);
    
    parfor (k=1:nc, par)
        tmpf = zeros(lat, 'single');
        tmpc = zeros(lat, 'single');
        for n=1:numel(f)
            tmpf = tmpf + single(f{n}(:,:,:,k));
            tmpc = tmpc + single(numeric(c{n}));
        end
        tmpf = tmpf ./ tmpc;
        mu(:,:,:,k) = tmpf;
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
    parfor (z=1:dim(3), par)
        tmpf = zeros([lat(1:2) 1 nc], 'single');
        for n=1:numel(f)
            tmpf = tmpf + single(f{n}(:,:,z,:));
        end
        mu(:,:,z,:) = tmpf ./ tmpc(:,:,z);
    end
    
    if ~isempty(output)
        mu = saveOnDisk(output, mu, 'name', 'mu');
    end
end

function mu = loopNone(f, c, output)

    dim = [size(f{1}) 1 1];
    lat = dim(1:3);
    nc  = dim(4);
    
    mu   = zeros([lat nc], 'single');
    tmpc = zeros(lat, 'single');
    for n=1:numel(f)
        mu   = mu + single(numeric(f{n}));
        tmpc = tmpc + single(numeric(c{n}));
    end
    mu = mu ./ tmpc;
    
    if ~isempty(output)
        mu = saveOnDisk(output, mu, 'name', 'mu');
    end
end