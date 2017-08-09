function saveArray(obj, s, ~)
% FORMAT obj.saveArray(s, ~)
% s - Name of the property which triggered the listener.
% 
% Event should be 'PreSet'

    p = s.Name;
    st = obj.disableListeners(p);
    if isa(obj.(p), 'nifti')
        try
            obj.private = obj.(p);
        catch e
            warning('Callback from %s: savearray failed\n', p);
            obj.enableListeners(st, p);
            rethrow(e)
        end
    end
    obj.enableListeners(st, p);
end