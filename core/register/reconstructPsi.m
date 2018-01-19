function psi = reconstructPsi(A, phi, varargin)
% FORMAT psi = reconstructPsi(A, phi, ('lat', lat), ('Mf', Mf), ('Mmu', Mmu))
%
% ** Required **
% A    - Affine transform (warps f to mu).
% phi  - Direct diffeomorphism (warps f to mu).
% ** Optional **
% lat  - Output lattice [default: same as phi]
% Mf   - Image voxel-to-world affine mapping.
% Mmu  - Template voxel-to-world affine mapping.
% ** Output **
% psi  - Complete direct transform: psi = A(phi)
%
% Reconstruct the inverse transform.

    % --- Parse inputs
    p = inputParser;
    p.FunctionName = 'reconstructPsi';
    p.addRequired('A',    @(X) issame(size(X), [4 4]));
    p.addRequired('phi',  @checkarray);
    p.addParameter('lat', [], @checkarray);
    p.addParameter('Mf',  eye(4), @(X) issame(size(X), [4 4]));
    p.addParameter('Mmu', eye(4), @(X) issame(size(X), [4 4]));
    p.addParameter('output', []);
    p.addParameter('debug',  false);
    p.parse(A, phi, varargin{:});
    lat = p.Result.lat;
    Mf  = p.Results.Mf;
    Mmu = p.Results.Mmu;
    
    if p.Results.debug, fprintf('* reconstructPsi\n'); end;
    
    % --- Load data
    A   = numeric(A);
    phi = numeric(phi);
    
    % --- Reconstruct
    if isempty(lat)
        psi = spm_warps('compose', Mmu \ A \ Mf, phi);
    else
        lat = [size(phi) 1];
        lat = lat(1:3);
        id  = spm_warps('identity', lat);
        psi = spm_warps('compose', Mmu \ A \ Mf, phi, id);
    end
    
    % --- Write on disk
    if ~isempty(p.Results.output)
        psi = saveOnDisk(p.Results.output, psi, 'name', 'psi');
    end
    
end