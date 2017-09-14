function gmu = templateGrad(mu, varargin)
% FORMAT gmu = templateGrad(mu, (itrp), (bnd), ...)
%
% ** Required **
% mu   - Template
% ** Optional **
% itrp - 1 or 3 Interpolation order [1]
% bnd  - 1 or 3 Boundary condition (0/1 = mirror/circulant) [1]
% ** Output **
% gmu  - Template spatial gradients
%
% Compute spatial gradients (w.r.t. local deformations) of the template
% image.

    % --- Parse inputs
    p = inputParser;
    p.FunctionName = 'templateGrad';
    p.addRequired('mu', @checkarray);
    p.addOptional('itrp', 1, @isnumeric);
    p.addOptional('bnd',  1, @isnumeric);
    p.addParameter('output', []);
    p.addParameter('debug', false);
    p.parse(mu, varargin{:});
    itrp   = p.Results.itrp;
    bnd    = p.Results.bnd;
    output = p.Results.output;
    debug  = p.Results.debug;
    
    
    if debug, fprintf('* templateGrad\n'); end;
    
    % --- Default parameters
    if numel(itrp) == 1
        itrp = itrp * [double(size(mu, 1) > 1) ...
                       double(size(mu, 2) > 1) ...
                       double(size(mu, 3) > 1)];
    end
    if numel(bnd) == 1
        bnd = [bnd bnd bnd];
    end

    % --- Read dimensions
    dim = [size(mu) 1 1 1 1];
    dim = dim(1:4);
    lat = dim(1:3);  % Template lattice dimension
    nc  = dim(4);    % Number of classes/modalities
    id  = single(warps('identity', lat)); % Identity transform

    % --- Reserve space for the complete gradient images
    if ~isempty(output)
        gmu = prepareOnDisk(output, [dim 3]);
    else
        gmu = zeros([dim 3], 'single');
    end

    % --- Do it
    % 1) Load one class
    % 2) Interpolation
    % 3) Get gradients
    % 4) Store data
    for i=1:nc
        mu1 = single(mu(:,:,:,i));
        c = spm_diffeo('bsplinc', mu1, [itrp bnd]);
        [~, G1, G2, G3] = spm_diffeo('bsplins', c, id, [itrp bnd]);
        gmu(:,:,:,i,1) = G1;
        gmu(:,:,:,i,2) = G2;
        gmu(:,:,:,i,3) = G3;
        clear c G1 G2 G3
    end

    % --- Write on disk
    if ~isempty(output)
        gmu = saveOnDisk(output, gmu, 'name', 'gmu');
    end
        
end