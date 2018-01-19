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
    
    if p.Results.debug, fprintf('* reconstructIPsi\n'); end;
    
    % --- Load data
    A    = numeric(A);
    iphi = numeric(iphi);
    
    if isempty(lat)
        lat = [size(iphi) 1];
        lat = lat(1:3);
    end
    
    % --- Reconstruct
    id   = spm_warps('identity', lat);
    ipsi = spm_warps('compose', iphi, Mmu \ (A \ Mf), id);
    
    % --- Write on disk
    if ~isempty(p.Results.output)
        ipsi = saveOnDisk(p.Results.output, ipsi, 'name', 'ipsi');
    end
    
end