function [lat, Mmu, vs] = autoTemplateDim(dat, vs)
%__________________________________________________________________________
% Compute template dimensions from input dmiensions and voxel-to-world
% mappings.
%--------------------------------------------------------------------------
% FORMAT [lat, Mmu, vs] = autoTemplateDim(dat, (vs))
%
% REQUIRED
% --------
% dat - Structure array with fields f (image) and Mf (voxel-to-world)
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

    if isempty(vs)
        vs = [0 0 0];
        for i=1:numel(dat)
            vs = vs + sqrt(sum(dat(i).Mf(1:3,1:3).^2));
        end
        vs = vs/numel(dat);
    end
    
    Mmu = eye(4);
    Mmu([1 6 11]) = vs;
    Mmu(1) = -Mmu(1);
    
    minpos = [inf inf inf]';
    maxpos = [-inf -inf -inf]';
    meantr = [0 0 0]';
    for i=1:numel(dat)
        corners = ones(4, 8);
        corners(1:3,1) = [1 1 1]';
        corners(1:3,2) = [1 1 size(dat(i).f, 3)]';
        corners(1:3,3) = [1 size(dat(i).f, 2) 1]';
        corners(1:3,4) = [size(dat(i).f, 1) 1 1]';
        corners(1:3,5) = [1 size(dat(i).f, 2) size(dat(i).f, 3)]';
        corners(1:3,6) = [size(dat(i).f, 1) 1 size(dat(i).f, 3)]';
        corners(1:3,7) = [size(dat(i).f, 1) size(dat(i).f, 2) 1]';
        corners(1:3,8) = [size(dat(i).f, 1) size(dat(i).f, 2) size(dat(i).f, 3)]';
        for j=1:8
            pushed = Mmu \ dat(i).Mf * corners(:,j);
            minpos = min(minpos, pushed(1:3));
            maxpos = max(maxpos, pushed(1:3));
        end
        meantr = meantr + dat(i).Mf(1:3,4);
    end
    meantr = meantr / numel(dat);
    Mmu(1:3,4) = meantr;
    lat = ceil(maxpos - minpos + 1);
    lat = lat(1:3)';
end