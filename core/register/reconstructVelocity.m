function v = reconstructVelocity(varargin)
% FORMAT (v) = obj.reconstructVelocity(('latent', z), ('subspace', W,),
%                                      ('residual' r), ('sigma', s))
% ** Keyword arguments **
% z    - Latent coordinates [do not use]
% W    - Principal subspace [do not use]
% r    - Residual velocity field [do not use]
% s    - Variance captured by the residual field [1]
% loop - How to split computation ('none', 'slice', 'component') [auto]
% par  - If true, parallelise processing [false]
% ** Output **
% v    - Output velocity
%
% Reconstruct the velocity field from its latent PC coordinates:
% > v = W * z + s * r

    % --- Parse inputs
    p = inputParser;
    p.FunctionName = 'reconstructVelocity';
    p.addParameter('latent',   [],    @checkarray);
    p.addParameter('subspace', [],    @checkarray);
    p.addParameter('residual', [],    @checkarray);
    p.addParameter('sigma',    1,     @isscalar);
    p.addParameter('loop',     '',    @(X) ischar(X) && any(strcmpi(X, {'slice', 'none', ''})));
    p.addParameter('par',      false, @isscalar);
    p.addParameter('output',   []);
    p.addParameter('debug',    false, @isscalar);
    p.parse(varargin{:});
    z      = p.Results.latent;
    W      = p.Results.subspace;
    r      = p.Results.residual;
    s      = p.Results.sigma;
    par    = p.Results.par;
    loop   = p.Results.loop;
    output = p.Results.output;
    debug  = p.Results.debug;
    
    if debug, fprintf('* reconstructVelocity\n'); end
    
    % --- Check that there are enough arguments
    if checkarray(z) + checkarray(W) == 1
        error('latent and subpsace keywords must be used together')
    end
    if ~checkarray(z) && ~checkarray(r)
        error('No argument provided')
    end
    
    % --- Dim info
    if checkarray(W)
        dim = [size(W) 1 1 1];
        lat = dim(1:3);
        nz  = dim(3);
        nv  = dim(4);
        nk  = dim(5);
        isfa = isa(W, 'file_array');
        if numel(z) ~= nk
            warning('Z and W have different number of classes')
        end
    else
        dim = [size(r) 1 1 1];
        lat = dim(1:3);
        nz  = dim(3);
        nv  = dim(4);
        nk = 1;
        isfa = isa(r, 'file_array');
    end
    
    % --- Optimise parallelisation and splitting schemes
    [par, loop] = autoParLoop(par, loop, isfa, nz, nk);
    
    % --- Allocate output
    v = prepareOnDisk(output, [lat nv]);
    v(:,:,:,:) = 0;
    
    % --- Compute linear combination
    if checkarray(z)
        if strcmpi(loop, 'none')
            v(:,:,:,:) = lat2vel(z, W);
        elseif strcmpi(loop, 'slice')
            if ~par
                for iz=1:nz
                    v(:,:,iz,:) = lat2vel(z, W(:,:,iz,:,:));
                end
            elseif isa(W, 'file_array')
                parfor (iz=1:nz, par)
                    v(:,:,iz,:) = wrap_lat2vel(z, W, iz, 3);
                end
            else
                parfor (iz=1:nz, par)
                    v(:,:,iz,:) = lat2vel(z, W(:,:,iz,:,:));
                end
            end
        else % strcmpi(loop, 'component')
            tmp = zeros(size(v), 'single');
            if ~par
                for k=1:nk
                    tmp = tmp + lat2vel(z(k), W(:,:,:,:,k));
                end
            elseif isa(W, 'file_array')
                parfor (k=1:nk, par)
                    tmp = tmp + wrap_lat2vel(z, W, k, 5);
                end
            else
                parfor (k=1:nk, par)
                    tmp = tmp + lat2vel(z(k), W(:,:,:,:,k));
                end
            end
            v(:,:,:,:) = tmp(:,:,:,:);
            clear tmp
        end
    end

    % --- Add residual field
    if checkarray(r)
        if any(strcmpi(loop, {'none', 'component'}))
            v(:) = v(:) + s * r(:);
        else % strcmpi(loop, 'slice')
            if ~par
                for iz=1:nz
                    v(:,:,iz,:) = v(:,:,iz,:) + s * r(:,:,iz,:);
                end
            elseif isa(v, 'file_array') && isa(r, 'file_array')
                parfor (iz=1:nz, par)
                    wrap_addtov(v, iz, 3, s * slicevol(r, iz, 3));
                end
            elseif isa(v, 'file_array')
                parfor (iz=1:nz, par)
                    wrap_addtov(v, iz, 3, s * r(:,:,iz,:));
                end
            elseif isa(r, 'file_array')
                parfor (iz=1:nz, par)
                    v(:,:,iz,:) = v(:,:,iz,:) + s * slicevol(r, iz, 3);
                end
            else
                parfor (iz=1:nz, par)
                    v(:,:,iz,:) = v(:,:,iz,:) + s * r(:,:,iz,:);
                end
            end
        end
    end

    
    % --- Write on disk
    v = saveOnDisk(output, v);
end

% PERFORM W*Z
% -----------

function v = lat2vel(z, W)
    dim = [size(W) 1 1 1];
    lat = dim(1:3);
    nv = dim(4);
    nk = dim(5);
    
    % Check if Z/W dimensions are the same
    if numel(z) > nk
        z = z(1:nk);
    elseif numel(z) < nk
        W = W(:,:,:,:,1:numel(z));
    end
    nk = numel(z);
    
    v = reshape(reshape(numeric(W), [], nk) * z(:), [lat nv]);
end

% WRAPPERS FOR EFFICIENT PARFOR
% -----------------------------

function v = wrap_addtov(v, ind, dim, value)
    if nargin < 3
        dim = 0;
    end
    switch dim
        case 0
        case 1
            v(ind,:,:,:) = v(ind,:,:,:) + value;
        case 2
            v(:,ind,:,:) = v(:,ind,:,:) + value;
        case 3
            v(:,:,ind,:) = v(:,:,ind,:) + value;
        otherwise
            error('Cannot split v along dimension %d', dim)
    end
end

function v = wrap_lat2vel(z, W, ind, dim)
    if nargin < 4
        dim = 0;
    end
    switch dim
        case 0
        case 1
            W = W(ind,:,:,:,:);
        case 2
            W = W(:,ind,:,:,:);
        case 3
            W = W(:,:,ind,:,:);
        case 5
            W = W(:,:,:,:,ind);
            z = z(ind);
        otherwise
            error('Cannot split W along dimension %d', dim)
    end
    v = lat2vel(z, W, ind, dim);
end