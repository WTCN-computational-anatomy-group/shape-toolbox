function s = struct(obj)
    obj.disableListeners();
    s = builtin('struct', obj);
    obj.enableListeners();
end