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
    
    if debug, fprintf('* reconstructVelocity\n'); end;
    
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
            parfor (iz=1:nz, par)
                v(:,:,iz,:) = v(:,:,iz,:) + lat2vel(z, W(:,:,iz,:,:));
            end
        else % strcmpi(loop, 'component')
            tmp = zeros(size(v), 'single');
            parfor (k=1:nk, par)
                tmp = tmp + lat2vel(z(k), W(:,:,:,:,k));
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
            parfor (iz=1:nz, par)
                v(:,:,iz,:) = v(:,:,iz,:) + s * r(:,:,iz,:);
            end
        end
    end

    
    % --- Write on disk
    v = saveOnDisk(output, v);
end

function v = lat2vel(z, W)
    dim = [size(W) 1 1 1];
    lat = dim(1:3);
    nv = dim(4);
    nk = dim(5);
    
    v = reshape(reshape(numeric(W), [], nk) * z(:), [lat nv]);
end