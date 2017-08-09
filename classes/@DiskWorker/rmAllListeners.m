function rmAllListeners(obj)
% FORMAT obj.rmAllListeners()
%
% Deletes all listeners. This function should not be used usually
% (prefer disableListeners and enableListeners to partially freeze
% one or several listeners)

    cellfun(@delete, obj.alllisteners);
    obj.alllisteners = {};
    fields = fieldnames(obj.listeners);
    for i=1:numel(fields)
        obj.listeners.(fields{i}) = {};
    end
    
end