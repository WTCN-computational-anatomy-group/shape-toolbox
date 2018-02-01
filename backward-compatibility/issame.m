function bool = issame(a, b)
    if exist('issame', 'builtin')
        bool = builtin('issame', a, b);
        return
    end
    bool = a == b;
    bool = all(bool(:));
end