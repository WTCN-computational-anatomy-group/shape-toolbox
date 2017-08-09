function ok = checkarray(obj, name)
% FORMAT ok = obj.checkarray(name)
% obj   - DiskWorker object
% name  - Private name of the array to check
% Checks that the array dimensions are not null, that the file exists and
% that the content is up-to-date
%
% FORMAT ok = obj.checkarray(array)
% obj   - DiskWorker object
% array - Array or file_array to check
% Checks that the array dimensions are not null and, if a file_array,
% that the file exists.

    if ischar(name)
        ok = obj.utd.(name) ...
             && prod(size(obj.(name))) ~= 0 ...
             && exist(obj.(name).fname, 'file');
    elseif isa(name, 'file_array')
        ok = prod(size(name)) ~= 0 && exist(name.fname, 'file');
    else
        ok = ~isempty(name);
    end
        
end

