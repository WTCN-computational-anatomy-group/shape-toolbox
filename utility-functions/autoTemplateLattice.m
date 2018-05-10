function [lat, Mmu, vs] = autoTemplateLattice(dims, M, vs)
%__________________________________________________________________________
% Compute template dimensions from input dimensions and voxel-to-world
% mappings.
%--------------------------------------------------------------------------
% FORMAT [lat, Mmu, vs] = autoTemplateLattice(dims, M, (vs))
%
% REQUIRED
% --------
% dims - List (cell array) of all input dimensions
% M    - List (cell array) of all input voxel-to-world mappings
%
% OPTIONAL
% --------
% vs  - Template voxel size [mean voxel size of inputs]
%
% OUTPUT
% ------
% lat - Template dimension that englobes all input images
% Mmu - Template voxel-to-world
% vs  - Template voxel size
%__________________________________________________________________________

    if nargin < 2
        vs = [];
    end

    % ---------------------------------------------------------------------
    % Check input
    if ~iscell(M)
        M = {M};
    end
    if numel(M) == 1
        M = repmat(M, 1, numel(dims));
    end
    if numel(dims) ~= numel(M)
        error('There should be as many voxel-to-world (%d) as dimensions(%d)', numel(M), numel(dims))
    end
    
    % ---------------------------------------------------------------------
    % Compute mean voxel size
    if isempty(vs)
        vs = [0 0 0];
        count = 0;
        for i=1:numel(M)
            if isempty(M{i})
                continue
            end
            count = count + 1;
            vs    = vs + sqrt(sum(M{i}(1:3,1:3).^2));
        end
        vs = vs/count;
    end
    
    % ---------------------------------------------------------------------
    % Initialise template voxel-to-world
    Mmu           = eye(4);
    Mmu([1 6 11]) = vs;
    Mmu(1)        = -Mmu(1);
    
    % ---------------------------------------------------------------------
    % Map all bounding boxes to template space
    minpos = [inf inf inf]';
    maxpos = [-inf -inf -inf]';
    meantr = [0 0 0]';
    count  = 0;
    for i=1:numel(dims)
        if isempty(dims{i}) || prod(dims{i}) == 0
            continue
        end
        count = count + 1;
        dim = [dims{i} 1 1 1];
        corners        = ones(4, 8);
        corners(1:3,1) = [1 1 1]';
        corners(1:3,2) = [1 1 dim(3)]';
        corners(1:3,3) = [1 dim(2) 1]';
        corners(1:3,4) = [dim(1) 1 1]';
        corners(1:3,5) = [1 dim(2) dim(3)]';
        corners(1:3,6) = [dim(1) 1 dim(3)]';
        corners(1:3,7) = [dim(1) dim(2) 1]';
        corners(1:3,8) = [dim(1) dim(2) dim(3)]';
        for j=1:8
            pushed = Mmu \ M{i} * corners(:,j);
            minpos = min(minpos, pushed(1:3));
            maxpos = max(maxpos, pushed(1:3));
        end
        meantr = meantr + M{i}(1:3,4);
    end
    meantr     = meantr / count;
    Mmu(1:3,4) = meantr;
    lat = ceil(maxpos - minpos + 1);
    lat = lat(1:3)';
end