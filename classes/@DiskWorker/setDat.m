function setDat(obj, private_name, value)
    public_name = obj.dat2nii.(private_name);
    s = obj.disableListeners(public_name);
    obj.(public_name).dat = value;
    obj.enableListeners(s, public_name);
end