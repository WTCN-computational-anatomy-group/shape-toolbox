function [wf, gf] = warp(phi, f, varargin)
% _________________________________________________________________________
%
% Generic util for applying a transform to an image
%
% -------------------------------------------------------------------------
% FORMAT [(wf), (gf)] = warp(phi, f, (itrp), (bnd), ...)
%
% REQUIRED
% --------
% phi  - Transform [nx ny nz 3]
% f    - Image     [mx my mz k]
%
% OPTIONAL
% --------
% itrp - Interpolation: 1 or 3 interpolation order [1].
% bnd  - Boundary: 1 or 3 boundary conditions [1].
% grad - Only compute warped gradients (not the warped image)
%
% KEYWORD ARGUMENTS
% -----------------
% par  - If true, parallelise processing [false]
%
% OUTPUT
% ------
% wf   - Warped image     [nx ny nz k]
% gf   - Warped gradients [nx ny nz k 3]
% _________________________________________________________________________

    % --- Parse inputs
    p = inputParser;
    p.FunctionName = 'warp';
    p.addRequired('phi', @checkarray);
    p.addRequired('f',   @checkarray);
    p.addOptional('itrp',      1,         @isnumeric);
    p.addOptional('bnd',       1,         @isnumeric);
    p.addOptional('grad',      false,     @islogical);
    p.addParameter('par',      false,     @isscalar);
    p.addParameter('output',   []);
    p.addParameter('debug',    false,     @isscalar);
    p.parse(phi, f, varargin{:});
    itrp   = p.Results.itrp;
    bnd    = p.Results.bnd;
    par    = p.Results.par;
    grad   = p.Results.grad;
    output = p.Results.output;
    debug  = p.Results.debug;
    
    if debug, fprintf('* warp\n'); end
    
    % --- Default parameters
    if numel(itrp) == 1
        itrp = itrp * [double(size(f, 1) > 1) ...
                       double(size(f, 2) > 1) ...
                       double(size(f, 3) > 1)];
    end
    if numel(bnd) == 1
        bnd = [bnd bnd bnd];
    end
    
    % --- Allocate output
    if ~iscell(output)
        output = {output};
    end
    if numel(output) < 2
        output = [output {[]}];
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
    compute_image = ~grad;
    compute_grad  = nargout > 1 || grad;
    if compute_image
        wf = prepareOnDisk(output{1}, [lat nc]);
    end
    if compute_grad
        if compute_image
            gf = prepareOnDisk(output{2}, [lat nc 3]);
        else
            gf = prepareOnDisk(output{1}, [lat nc 3]);
        end
    end
    
    % --- Interpolate
    if compute_grad
        if compute_image
            for k=1:nc
                c = spm_diffeo('bsplinc', single(f(:,:,:,k)), [itrp bnd]);
                    [wf(:,:,:,k), gf(:,:,:,k,1), gf(:,:,:,k,2), gf(:,:,:,k,3)] ...
                        = spm_diffeo('bsplins', c, single(numeric(phi)), [itrp bnd]);
            end
        else
            for k=1:nc
                c = spm_diffeo('bsplinc', single(f(:,:,:,k)), [itrp bnd]);
                    [~, gf(:,:,:,k,1), gf(:,:,:,k,2), gf(:,:,:,k,3)] ...
                        = spm_diffeo('bsplins', c, single(numeric(phi)), [itrp bnd]);
            end
        end
    elseif ~par
        for k=1:nc
            c = spm_diffeo('bsplinc', single(f(:,:,:,k)), [itrp bnd]);
            wf(:,:,:,k) = spm_diffeo('bsplins', c, single(numeric(phi)), [itrp bnd]);
        end
    elseif isa(f, 'file_array')
        parfor (k=1:nc, par)
            c = spm_diffeo('bsplinc', single(slicevol(f, k, 4)), [itrp bnd]);
            wf(:,:,:,k) = spm_diffeo('bsplins', c, single(numeric(phi)), [itrp bnd]);
        end
    else
        parfor (k=1:nc, par)
            c = spm_diffeo('bsplinc', single(f(:,:,:,k)), [itrp bnd]);
            wf(:,:,:,k) = spm_diffeo('bsplins', c, single(numeric(phi)), [itrp bnd]);
        end
    end
    
    % --- Write on disk
    if compute_image
        wf = saveOnDisk(output{1}, wf);
    end
    if compute_grad 
        if compute_image
            gf = saveOnDisk(output{2}, gf);
        else
            gf = saveOnDisk(output{1}, gf);
        end
    end
    if grad
        wf = [];
        [gf, wf] = deal(wf, gf);
    end
end