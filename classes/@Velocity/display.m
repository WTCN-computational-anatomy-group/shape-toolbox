function display(obj)
    obj.disableListeners();
    builtin('disp', obj);
    obj.enableListeners();
end