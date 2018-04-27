function val = defval(a, fields, default)

    try
        val = eval(['a' fields]);
    catch
        val = default;
    end
        

end