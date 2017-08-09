function value = getDat(obj, private_name)
    public_name = obj.dat2nii.(private_name);
    s = obj.disableListeners(public_name);
    value = obj.(public_name).dat;
    obj.enableListeners(s, public_name);
end