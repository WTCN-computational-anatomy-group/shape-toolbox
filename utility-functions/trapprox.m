function t = trapprox(varargin)
% _________________________________________________________________________
%
%      Stochastic approximation of tr(A^m), where A=(H+L)\L or A=H+L.
%
% -------------------------------------------------------------------------
% FORMAT t = trapprox(L, ...)
% FORMAT t = trapprox(prm, ...)
% FORMAT t = trapprox(prm, H, ...)
%
% MANDATORY
% --------
% L       - Spare or full matrix
%    or
% prm     - Parameters of the prior precision (if spm_diffeo/spm_field)
%
% OPTIONAL
% --------
% H       - Hessian of the data term:         nx*ny*nz*6 or [empty/0]
%
% KEYWORD ARGUMENTS
% -----------------
% vs      - Lattice voxel size                              [1 1 1]
% fmg     - Full Multigrid options                          [cyc rlx]
%           cyc - Number of Full Multigrid cycles           [2]
%           rlx - Number of relaxation iterations per cycle [2]
% moments - Number of moments                               [1]
% samples - Number of random samples                        [10]
% matrix  - Matrix form                               'H+L'/['(H+L)\L']
% method  - Sampling method                      'gaussian'/['rademacher']
% type    - Input type                              'field'/['diffeo']
%
% OUTPUT
% ------
% t       - List of moments
%
% -------------------------------------------------------------------------
%
% It should be more precise than spm_diffeo('trapprox') but takes a bit
% longer. Speedup can be obtained by reducing the number of samples:
% empirically, as few as 10 samples yield a better approximation than
% spm_diffeo.
% _________________________________________________________________________


    % ---------------------------------------------------------------------
    % Parse inputs
    % ---------------------------------------------------------------------
    p = inputParser;
    p.addRequired('L');
    p.addOptional('H',        []);
    p.addParameter('vs',      [1 1 1]);
    p.addParameter('fmg',     [2 2]);
    p.addParameter('moments', 1);
    p.addParameter('samples', 5);
    p.addParameter('matrix',  '(H+L)\L');
    p.addParameter('method',  'rademacher');
    p.addParameter('type',    'diffeo');
    p.addParameter('dim',     []);
    p.parse(varargin{:});
    
    L     = p.Results.L;
    H     = p.Results.H;
    fmg   = p.Results.fmg;
    nm    = p.Results.moments;
    ns    = p.Results.samples;
    dim   = p.Results.dim;
    which = strcmpi(p.Results.matrix, 'H+L');
    samp  = str2func(p.Results.method);
    if strcmpi(p.Results.type, 'field')
        spm_func = @spm_field;
    else
        spm_func = @spm_diffeo;
    end

    % matrix or paramter form?
    if size(L,1) ~= size(L,2)
        prm = L;
        L   = [];
        prm = [p.Results.vs prm];
    end
    
    % empty Hessian?
    if isempty(H) || numel(H) == 1
        if isempty(dim)
            dim = [size(L, 1) 1];
        end
        which = true; % 'H+L' case
        no_hessian = true;
    else
        dim = [size(H) 1 1];
        lat = dim(1:3);
        [~,nc]  = spm_matcomp('SymIndices', dim(4), 'k');
        dim = [lat nc];
        H = single(H);
        no_hessian = false;
    end

    
    % ---------------------------------------------------------------------
    % Compute moments
    % ---------------------------------------------------------------------
    t = zeros(nm,1);
    for i=1:ns                      % samples
        v = samp(dim);
        m = v;
        for j=1:nm                  % moments
            % Compute A^m * v
            if which                % H+L
                if isempty(L)
                    if no_hessian
                        m = spm_func('vel2mom', single(m), double(prm));
                    else
                        m = spm_func('vel2mom', single(m), double(prm)) ...
                            + spm_matcomp('pointwise',H,m,true);
                    end
                else
                    if isempty(H) || numel(H) == 1
                        m = reshape(single(L*double(m(:))), size(m));
                    else
                        m = reshape(single(L*double(m(:))), size(m)) ...
                            + spm_matcomp('pointwise',H,m,true);
                    end
                end
            else                    % (H+L)\L
                if isempty(L)
                    m = spm_func('vel2mom', single(m), double(prm));
                else
                    m = reshape(single(L*double(m(:))), size(m));
                end
                m = spm_func('fmg', single(H), single(m), double([prm fmg]));
            end
            % Compute v' * (A^m * v)
            t(j) = t(j) + double(v(:)' * m(:));
        end
    end
    t = t/ns;

end

% -------------------------------------------------------------------------
% Random functions
% -------------------------------------------------------------------------
    
function x = rademacher(dim)
    x = single(binornd(1, 0.5, dim));
    x(x==0) = -1;
end

function x = gaussian(dim)
    x = single(normrnd(0, 1, dim));
end