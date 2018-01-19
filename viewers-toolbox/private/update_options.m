function optout = update_options(optin, new_opt)
    optout = optin;
    up_fnames = fieldnames(new_opt);
    for i=1:length(up_fnames)
        optout.(up_fnames{i}) = new_opt.(up_fnames{i});
    end
end