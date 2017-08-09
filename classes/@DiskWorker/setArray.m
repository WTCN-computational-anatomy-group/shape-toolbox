function setArray(obj, s, ~)
% FORMAT obj.setArray(s, ~)
% s - Name of the property which triggered the listener.
% 
% Event should be 'PostSet'

    p = s.Name;
    st = obj.disableListeners(p);
    try
        obj.(p) = setNifti(obj.(p), obj.private);
        obj.statusChanged(obj.nii2dat.(p));
        obj.private = [];
    catch e
        warning('Callback from %s: setarray failed\n', p);
        obj.enableListeners(st, p);
        rethrow(e)
    end
    obj.enableListeners(st, p);
end