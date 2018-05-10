function gmu = templateGrad(mu, varargin)
% _________________________________________________________________________
%
% Compute spatial gradients (w.r.t. local deformations) of the template
% image.
%
% -------------------------------------------------------------------------
% FORMAT gmu = templateGrad(mu, (bnd), ...)
%
% REQUIRED
% --------
% mu   - (log)-Template
%
% OPTIONAL
% --------
% bnd    - 1 or 3 Boundary condition (0/1 = mirror/circulant)   [1]
% output - file_array where to store the output                 []
% debug  - Debugging talk                                       [false]
%
% OUTPUT
% ------
% gmu  - Template spatial gradients
%
% -------------------------------------------------------------------------
% Gradients are computed at the control points. In this context, the 1st
% order is ill-suited (even if it is the order used when actually
% interpolating the template). Here, we use second order interpolation in
% the gradient direction and zero-th order in the other directions. This
% yields gradients that are coherent with those that would be returned by
% Matlab's `gradient` function.
% _________________________________________________________________________

    % --- Parse inputs
    p = inputParser;
    p.FunctionName = 'templateGrad';
    p.addRequired('mu', @checkarray);
    p.addOptional('bnd',  1, @(X) isnumeric(X) || islogical(X));
    p.addParameter('output', []);
    p.addParameter('debug', false);
    p.parse(mu, varargin{:});
    bnd    = p.Results.bnd;
    output = p.Results.output;
    debug  = p.Results.debug;
    
    
    if debug, fprintf('* templateGrad\n'); end
    
    % --- Default parameters
    if numel(bnd) == 1
        bnd = [bnd bnd bnd];
    end

    % --- Read dimensions
    dim = [size(mu) 1 1 1 1];
    dim = dim(1:4);
    lat = dim(1:3);  % Template lattice dimension
    nc  = dim(4);    % Number of classes/modalities
    id  = single(spm_warps('identity', lat)); % Identity transform

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
        [~,gmu(:,:,:,i,1),~,~] = spm_diffeo('bsplins', mu1, id, [2 0 0 bnd]);
        [~,~,gmu(:,:,:,i,2),~] = spm_diffeo('bsplins', mu1, id, [0 2 0 bnd]);
        [~,~,~,gmu(:,:,:,i,3)] = spm_diffeo('bsplins', mu1, id, [0 0 2 bnd]);
    end

    % --- Write on disk
    if ~isempty(output)
        gmu = saveOnDisk(output, gmu, 'name', 'gmu');
    end
        
end