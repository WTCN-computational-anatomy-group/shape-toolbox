function wf = warp(phi, f, varargin)
% FORMAT wf = warp(phi, f, (itrp), (bnd), ('par', par))
%
% ** Required **
% phi  - Transform
% f    - Image
% ** Optional
% itrp - Interpolation: 1 or 3 interpolation order [1].
% bnd  - Boundary: 1 or 3 boundary conditions [1].
% ** Keyword arguments **
% par  - If true, parallelise processing [false]
% ** Output **
% wf   - Warped image
%
% Generic util for applying a transform to an image

    % --- Parse inputs
    p = inputParser;
    p.FunctionName = 'warp';
    p.addRequired('phi', @checkarray);
    p.addRequired('f',   @checkarray);
    p.addOptional('itrp',      1, @isnumeric);
    p.addOptional('bnd',       1, @isnumeric);
    p.addParameter('par',      false,     @isscalar);
    p.addParameter('output',   []);
    p.addParameter('debug',    false,     @isscalar);
    p.parse(phi, f, varargin{:});
    itrp   = p.Results.itrp;
    bnd    = p.Results.bnd;
    par    = p.Results.par;
    output = p.Results.output;
    debug  = p.Results.debug;
    
    if debug, fprintf('* warp\n'); end;
    
    % --- Default parameters
    if numel(itrp) == 1
        itrp = itrp * [double(size(f, 1) > 1) ...
                       double(size(f, 2) > 1) ...
                       double(size(f, 3) > 1)];
    end
    if numel(bnd) == 1
        bnd = [bnd bnd bnd];
    end
    
    % --- Optimise parallelisation and splitting schemes
    par = autoParLoop(par, 'component', isa(f, 'file_array'), 1, size(f,4));
    
    % --- Dim info
    dim = [size(f) 1 1];
    dim = dim(1:4);
    nc  = dim(4);               % Number of classes/modalities
    dim = [size(phi) 1 1];
    dim = dim(1:4);
    lat = dim(1:3);             % Dimension of the output lattice
    
    % --- Allocate output
    wf = prepareOnDisk(output, [lat nc]);
    
    % --- Interpolate
    parfor (k=1:nc, par)
        c = spm_diffeo('bsplinc', single(f(:,:,:,k)), [itrp bnd]);
        wf(:,:,:,k) = spm_diffeo('bsplins', c, single(numeric(phi)), [itrp bnd]);
    end
    
    % --- Write on disk
    wf = saveOnDisk(output, wf);
end