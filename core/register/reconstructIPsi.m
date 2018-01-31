function ipsi = reconstructIPsi(A, iphi, varargin)
% FORMAT ipsi = reconstructIPsi(A, iphi, ('lat', lat), ('Mf', Mf), ('Mmu', Mmu))
%
% ** Required **
% A    - Affine transform (warps f to mu).
% iphi - Inverse diffeomorphism (warps mu to f).
% ** Optional **
% lat  - Output lattice [default: same as iphi]
% Mf   - Image voxel-to-world affine mapping.
% Mmu  - Template voxel-to-world affine mapping.
% ** Output **
% ipsi - Complete inverse transform: ipsi = iphi(inv(A))
%
% Reconstruct the inverse transform.

    % --- Parse inputs
    p = inputParser;
    p.FunctionName = 'reconstructIPsi';
    p.addRequired('A',    @(X) issame(size(X), [4 4]));
    p.addRequired('iphi', @checkarray);
    p.addParameter('lat', [], @checkarray);
    p.addParameter('Mf',  eye(4), @(X) issame(size(X), [4 4]));
    p.addParameter('Mmu', eye(4), @(X) issame(size(X), [4 4]));
    p.addParameter('output', []);
    p.addParameter('debug',  false);
    p.parse(A, iphi, varargin{:});
    lat = p.Results.lat;
    Mf  = p.Results.Mf;
    Mmu = p.Results.Mmu;
    
    if p.Results.debug, fprintf('* reconstructIPsi\n'); end
    
    % --- Load data
    A    = numeric(A);
    iphi = numeric(iphi);
    
    if isempty(lat)
        lat = size(iphi);
    end
    lat = [lat 1 1 1];
    lat = lat(1:3);
    
    % --- Reconstruct
    % Assume phi is encoded by splines coeffs
    ixi  = spm_warps('affine', Mmu \ (A \ Mf), [lat 3]);
    id   = single(spm_warps('identity', size(iphi)));
    ipsi = single(spm_warps('identity', size(ixi)));
    for d=1:size(ipsi, 4)
        ipsi(:,:,:,d) = ixi(:,:,:,d) + spm_diffeo('bsplins', ...
            single(iphi(:,:,:,d)-id(:,:,:,d)), single(ixi), [1 1 1  1 1 1]);
    end
    clear ixi id iphi
    
    % --- Write on disk
    if ~isempty(p.Results.output)
        ipsi = saveOnDisk(p.Results.output, ipsi, 'name', 'ipsi');
    end
    
end