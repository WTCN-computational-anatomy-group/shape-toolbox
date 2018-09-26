function val = defval(a, fields, default)
% FORMAT val = defval(a, fields, default)
% a       - structure
% fields  - string representing the (hierarchical) fields to access
% default - Value to return is the fields do not exist.
%
% EXAMPLE status = defval(result, 'output.status', 0);
% Will return > result.output.status    if it exists
%             > 0                       if it does not exist

    if fields(1) ~= '.'
        fields = ['.' fields];
    end

    try
        val = eval(['a' fields]);
    catch
        val = default;
    end
        

end