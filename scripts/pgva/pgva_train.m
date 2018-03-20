function pgva_train(jsn_input,jsn_opt)
% FORMAT pgva_train(jsn_input,jsn_opt)
% jsn_input - Json file containing input filenames
% jsn_opt   - Json file containing option values
%
% Train a PGVA shape model on some data. This is a wrapper around the
% `pgva_model` function that takes JSON files as inputs and that can be
% compiled and used on command line.
%
% Both json should contain structured values (equivalent to a a maltab
% structure). Type `help pgva_model` for details on these structures and
% their fields.

    if nargin > 1
        opt = spm_jsonread(jsn_opt);
    else
        opt = struct;
    end
    
    if nargin == 0 || (ischar(jsn_input) && ...
            (strcmp(jsn_input,'--help') || ...
             strcmp(jsn_input,'--h')    || ...
             strcmp(jsn_input,'-h')))
        show_instructions;
    else
        input = spm_jsonread(jsn_input);
        if isfield(opt, 'z') && isfield(opt.z, 'A0')
            if ischar(opt.z.A0) && exist(opt.z.A0, 'file')
                opt.z.A0 = load(opt.z.A0);
            elseif isnumeric(opt.z.A0) && isvector(opt.z.A0)
                opt.z.A0 = diag(opt.z.A0);
            end
        end
        if isfield(opt, 'q') && isfield(opt.z, 'A0')
            if ischar(opt.q.A0) && exist(opt.q.A0, 'file')
                opt.q.A0 = load(opt.q.A0);
            elseif isnumeric(opt.q.A0) && isvector(opt.q.A0)
                opt.q.A0 = diag(opt.q.A0);
            end
        end
        pgva_model(input,opt);
    end

end

% Temporary help
function show_instructions
    help_str = '';
    disp('Usage: pgva_train input.jsn opt.jsn');
    disp(' ');
    disp(help_str);
    disp(' ');
end