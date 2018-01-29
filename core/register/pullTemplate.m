function [wa, ga] = pullTemplate(phi, a, varargin)
% _________________________________________________________________________
%
% Generic util for "pulling" a template into image space.
%
% -------------------------------------------------------------------------
% FORMAT [(wa), (ga)] = pullTemplate(ipsi, a, ...)
%
% REQUIRED
% --------
% ipsi  - Transform [nx ny nz 3]
% a     - Template  [mx my mz k]
%
% KEYWORD ARGUMENTS
% -----------------
% grad   - Only compute warped gradients (not the warped image) [false]
% par    - If true, parallelise processing over classes/modalities [false]
% output - Cell of file_array to use as output [return numeric array]
% debug  - Print debugging info [false]
%
% OUTPUT
% ------
% wa   - Warped template  [nx ny nz k]
% ga   - Warped gradients [nx ny nz k 3]
%
% -------------------------------------------------------------------------
% - Template are encoded as spline coefficients, so there is no need to
%   compute those first, contrary to what is done in `warp`.
% - Since the `pushImage` operation is only implemented for trilinear
%   inteprolation, we also restrict `pull` to 1st order.
% - Boundary conditions are circulant.
% _________________________________________________________________________

    % --- Parse inputs
    p = inputParser;
    p.FunctionName = 'warp';
    p.addRequired('phi', @checkarray);
    p.addRequired('a',   @checkarray);
    p.addParameter('grad',     false,     @islogical);
    p.addParameter('par',      false,     @isscalar);
    p.addParameter('output',   []);
    p.addParameter('debug',    false,     @isscalar);
    p.parse(phi, a, varargin{:});
    grad   = p.Results.grad;
    par    = p.Results.par;
    output = p.Results.output;
    debug  = p.Results.debug;
    
    if debug, fprintf('* pullTemplate\n'); end
    
    % --- Allocate output
    if ~iscell(output)
        output = {output};
    end
    if numel(output) < 2
        output = [output {[]}];
    end
    
    % --- Optimise parallelisation and splitting schemes
    par = autoParLoop(par, 'component', isa(a, 'file_array'), 1, size(a,4));
    
    % --- Dim info
    dim = [size(a) 1 1];
    dim = dim(1:4);
    nc  = dim(4);               % Number of classes/modalities
    dim = [size(phi) 1 1];
    dim = dim(1:4);
    lat = dim(1:3);             % Dimension of the output lattice
    
    % --- Allocate output
    compute_image = ~grad;
    compute_grad  = nargout > 1 || grad;
    if compute_image
        wa = prepareOnDisk(output{1}, [lat nc]);
    end
    if compute_grad
        if compute_image
            ga = prepareOnDisk(output{2}, [lat nc 3]);
        else
            ga = prepareOnDisk(output{1}, [lat nc 3]);
        end
    end
    
    % --- Interpolate
    if compute_grad
        if compute_image
            for k=1:nc
                [wa(:,:,:,k), ga(:,:,:,k,1), ga(:,:,:,k,2), ga(:,:,:,k,3)] ...
                    = spm_diffeo('bsplins', single(a(:,:,:,k)), single(numeric(phi)), [itrp bnd]);
            end
        else
            for k=1:nc
                [~, ga(:,:,:,k,1), ga(:,:,:,k,2), ga(:,:,:,k,3)] ...
                    = spm_diffeo('bsplins', single(a(:,:,:,k)), single(numeric(phi)), [itrp bnd]);
            end
        end
    elseif ~par
        for k=1:nc
            wa(:,:,:,k) = spm_diffeo('bsplins', single(a(:,:,:,k)), single(numeric(phi)), [itrp bnd]);
        end
    elseif isa(f, 'file_array')
        parfor (k=1:nc, par)
            wa(:,:,:,k) = spm_diffeo('bsplins', single(slicevol(a, k, 4)), single(numeric(phi)), [itrp bnd]);
        end
    else
        parfor (k=1:nc, par)
            wa(:,:,:,k) = spm_diffeo('bsplins', single(a(:,:,:,k)), single(numeric(phi)), [itrp bnd]);
        end
    end
    
    % --- Write on disk
    if compute_image
        wa = saveOnDisk(output{1}, wa);
    end
    if compute_grad 
        if compute_image
            ga = saveOnDisk(output{2}, ga);
        else
            ga = saveOnDisk(output{1}, ga);
        end
    end
    if grad
        wa = [];
        [ga, wa] = deal(wa, ga);
    end
end