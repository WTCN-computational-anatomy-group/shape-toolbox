function pgra_train(jsn_input,jsn_opt)
% FORMAT pgra_train(jsn_input,jsn_opt)
% jsn_input - Json file containing input filenames
% jsn_opt   - Json file containing option values
%
% Train a PGRA shape model on some data. This is a wrapper around the
% `pgra_model` function that takes JSON files as inputs and that can be
% compiled and used on command line.
%
% Both json should contain structured values (equivalent to a a maltab
% structure). Type `help pgra_model` for details on these structures and
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
        pgra_model(input,opt);
    end

end

% Temporary help
function show_instructions
    help_str = '';
    disp('Usage: pgra_train input.jsn opt.jsn');
    disp(' ');
    disp(help_str);
    disp(' ');
end