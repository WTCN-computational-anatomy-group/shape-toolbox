function b = bytes(obj)
    warning off 'MATLAB:structOnObject' 
    s = builtin('struct', obj); 
    w = whos('s'); 
    w.bytes
    warning on 'MATLAB:structOnObject'
end