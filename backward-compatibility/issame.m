function bool = issame(a, b)
    if exist('issame', 'builtin')
        bool = issame(a, b);
        return
    end
    bool = a == b;
    bool = all(bool(:));
end