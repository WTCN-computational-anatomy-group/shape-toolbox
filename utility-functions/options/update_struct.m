function out = update_struct(new, old)
% FORMAT out = update_struct(new, old)
% new - Struct (can be hierarchical) with new values.
% old - Struct (can be hierarchical) with original values.
%
% Updates the values of old with those of new.
% - Fields of old wchich are not in new are kept unchanged.
% - If some fields of new are also structures (hierarchical object), they 
%   are also updated and not entirely overwritten.
    if nargin < 1
        out = new;
    else
        out = old;
        newfn = fieldnames(new);
        for i=1:length(newfn)
            fn = newfn{i};
            if isstruct(new.(fn))
                if ~isfield(out, fn) || ~isstruct(out.(fn))
                    out.(fn) = struct;
                end
                out.(fn) = update_struct(new.(fn), out.(fn));
            else
                out.(fn) = new.(fn);
            end
        end
    end
end