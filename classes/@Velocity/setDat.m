function setDat(obj, private_name, value)
    if isfield(obj.dat2nii, private_name)
        public_name = obj.dat2nii.(private_name);
    else
        public_name = private_name;
    end
    s = obj.disableListeners(public_name);
    obj.(public_name).dat = value;
    obj.enableListeners(s, public_name);
end