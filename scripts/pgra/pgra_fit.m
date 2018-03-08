function pgra_fit(jsn_input,jsn_opt)
% FORMAT pgra_fit(jsn_input,jsn_opt)
% jsn_input - Json file containing input filenames
% jsn_opt   - Json file containing option values
%
% Fit a PGRA shape model on new data. This is a wrapper around the
% `pgra_model` function that takes JSON files as inputs and that can be
% compiled and used on command line.
%
% Both json should contain structured values (equivalent to a a maltab
% structure). Type `help pgra_model` for details on these structures and
% their fields.
%
% A few pgra_model options are overriden so that the model is in 'fit' 
% mode. In particular, optimise.pg, optimise.tpl, optimise.z.A, 
% optimise.r.l and optimise.q.A are all set to false (i.e., they are 
% fixed model parameters).

    if nargin > 1
        opt = spm_jsonread(jsn_opt);
    else
        opt = struct;
    end
    
    % Force values to forbid learning population parameters
    opt.optimise.pg  = false;
    opt.optimise.tpl = false;
    if isfield(opt.optimise, 'z')
        if isstruct(opt.optimise.z)
            opt.optimise.z.A = false;
        else
            value = opt.optimise.z;
            opt.optimise.z = struct;
            opt.optimise.z.A = false;
            opt.optimise.z.z = value;
        end
    end
    if isfield(opt.optimise, 'q')
        if isstruct(opt.optimise.q)
            opt.optimise.q.A = false;
        else
            value = opt.optimise.q;
            opt.optimise.q = struct;
            opt.optimise.q.A = false;
            opt.optimise.q.q = value;
        end
    end
    if isfield(opt.optimise, 'r')
        if isstruct(opt.optimise.r)
            opt.optimise.r.l = false;
        else
            value = opt.optimise.r;
            opt.optimise.r = struct;
            opt.optimise.r.l = false;
            opt.optimise.r.r = value;
        end
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
    disp('Usage: pgra_fit input.jsn opt.jsn');
    disp(' ');
    disp(help_str);
    disp(' ');
end