function getProperty(obj, s, ~)
% FORMAT obj.getProperty(s, ~)
% s - Name of the property which triggered the listener.
% 
% Event should be 'PreGet'
% Call the appropriate updater. If a nifti, call create (i.e., create file
% on disk if needed and write the header).

    p = s.Name;
    st = obj.disableListeners(p);
    try
        if isa(obj.(p), 'nifti') && isfield(obj.updaters, obj.nii2dat.(p))
            field = obj.nii2dat.(p);
            if isfield(obj.utd, field) ...
                    && ~obj.utd.(field) ...
                    && isfield(obj.updaters, field)
                obj.updaters.(field)(obj);
            end 
            create(obj.(p));
        elseif isfield(obj.utd, p) ...
                && ~obj.utd.(p) ...
                && isfield(obj.updaters, p)
            obj.updaters.(p)(obj);
        end 
    catch e
        warning('Callback from %s: getProperty failed\n', p);
        obj.enableListeners(st, p);
        rethrow(e)
    end
    obj.enableListeners(st, p);
end