function gmu = templateGrad(mu, varargin)
% FORMAT gmu = obj.computeTemplateGrad(mu, (itrp), (bnd))
% mu   - Template
% itrp - Interpolation order [default: 1 1 1]
% bnd  - Boundary condition (0/1 = mirror/circulant) [default: 1 1 1]
% gmu  - Template spatial gradients
%
% Compute spatial gradients (w.r.t. local deformations) of the template
% image.

    % --- Parse inputs
    p = inputParser;
    p.FunctionName = 'templateGrad';
    p.addRequired('mu', @checkarray);
    p.addOptional('itrp', [1 1 1]);
    p.addOptional('bnd', [1 1 1]);
    p.addParameter('output', []);
    p.addParameter('debug', false);
    p.parse(mu, varargin{:});
    
    if p.Results.debug, fprintf('* templateGrad\n'); end;

    % --- Read dimensions
    dim = [size(mu) 1 1 1 1];
    dim = dim(1:4);
    lat = dim(1:3);  % Template lattice dimension
    nc  = dim(4);    % Number of classes/modalities
    id  = single(warps('identity', lat)); % Identity transform

    % --- Reserve space for the complete gradient images
    if ~isempty(p.Results.output)
        gmu = prepareOnDisk(p.Results.output, [dim 3]);
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
        c = spm_diffeo('bsplinc', mu1, [p.Results.itrp p.Results.bnd]);
        [~, G1, G2, G3] = spm_diffeo('bsplins', c, id, [p.Results.itrp p.Results.bnd]);
        gmu(:,:,:,i,1) = G1;
        gmu(:,:,:,i,2) = G2;
        gmu(:,:,:,i,3) = G3;
        clear c G1 G2 G3
    end

    % --- Write on disk
    if ~isempty(p.Results.output)
        gmu = saveOnDisk(p.Results.output, gmu, 'name', 'gmu');
    end
        
end