function disp(obj)
    obj.disableListeners();
    builtin('disp', obj);
    obj.enableListeners();
end