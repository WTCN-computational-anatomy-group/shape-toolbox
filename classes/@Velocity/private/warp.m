function wf = warp(phi, f, bs, wf)
% FORMAT wf = warp(phi, f, bs, wf)
% phi - Transform
% f   - Image
% bs  - Interpolation options: 3 orders and 3 boundary conditions.
%       (See: help spm_bsplinc)
% wf  - Warped image
%
% Generic util for applying a transform to an image

    if nargin < 3 || isempty(bs)
        bs = [1 1 1 1 1 1]; % Linear/Circulant
    end
    
    % --- Dim info
    dim = [size(f) 1 1 1 1];
    dim = dim(1:4);
    nc  = dim(4);               % Number of classes/modalities
    dim = [size(phi) 1 1 1 1];
    dim = dim(1:4);
    lat = dim(1:3);             % Dimension of the output lattice
    
    % --- Prepare output
    if nargin < 4
        wf = zeros([lat nc], 'single');
    elseif isa(wf, 'file_array')
        wf.dim = [lat nc];
    end
    
    % --- Interpolate
    for k=1:nc
        % If file_array, load one class (saves memory)
        c = spm_diffeo('bsplinc', single(f(:,:,:,k)), bs);
        % Write directly to disk
        wf(:,:,:,k) = spm_diffeo('bsplins', c, single(numeric(phi)), bs);
    end
end