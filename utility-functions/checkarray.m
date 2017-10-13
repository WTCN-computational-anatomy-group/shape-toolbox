function ok = checkarray(array)
% FORMAT ok = checkarray(array)
% array - Array or file_array to check
%
% Checks that the array dimensions are not null and, if a file_array,
% that the file exists.

    if isa(array, 'file_array')
        ok = prod(size(array)) ~= 0 ...
            && all(cellfun(@(X) exist(X, 'file'), {array.fname}));
    elseif isnumeric(array)
        ok = ~isempty(array);
    else
        ok = false;
    end
        
end

