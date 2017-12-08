function a = slicevol(a, ind, dim)
% FORMAT s = slicevol(a, ind, dim)
% a   - array (numeric or file_array)
% ind - Index or range of indices
% dim - Dimension along which to slice
% s   - Selected slice(s)
%
% 
% FORMAT s = slicevol(a, {ind1, ind2, ind2, ...})
% a    - array (numeric or file_array)
% indx - Index or range of indices for dimension x
% s   - Selected subvolume
% 
% This function allows to select a subvolume from a volume. It is useful
% when the processing of a file_array is split and parallelised with a
% parfor.

    % bounding box
    if iscell(ind)
        if numel(ind) < 1
            ind{1} = 1:size(a, 1);
        end
        if numel(ind) < 2
            ind{2} = 1:size(a, 2);
        end
        if numel(ind) < 3
            ind{3} = 1:size(a, 3);
        end
        if numel(ind) < 4
            ind{4} = 1:size(a, 4);
        end
        if numel(ind) < 5
            ind{5} = 1:size(a, 5);
        end
        a = a(ind{1}, ind{2}, ind{3}, ind{4}, ind{5});
        return
    end
    % slice
    if nargin < 3
        dim = 0;
    end
    switch dim
        case 0
        case 1
            a = a(ind,:,:,:,:);
        case 2
            a = a(:,ind,:,:,:);
        case 3
            a = a(:,:,ind,:,:);
        case 4
            a = a(:,:,:,ind,:);
        case 5
            a = a(:,:,:,:,ind);
        otherwise
            error('Cannot slice along dimension %d', dim);
    end
end