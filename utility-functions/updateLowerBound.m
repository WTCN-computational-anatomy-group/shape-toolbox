function model = updateLowerBound(model, gain)
% FORMAT model = updateLowerBound(model, ('gain'))
% model - Model structure
% gain  - Last call before new EM loop ? [false]
%
% Compute the lower bound from its subparts
% + Compute gain against the previous VEM iteration if [last == true]
% + Some stuff to have nice plots where x values are actual EM iterations.

    if nargin < 2
        gain = '';
    end
    gain = strcmpi(gain, 'gain');
    
    if gain
        % Compute gain
        % ------------
        if ~isfield(model.lb.lb, 'gainlist')
            model.lb.lb.gainlist = [];
        end
        if isempty(model.lb.lb.list)
            model.lb.lb.gain = inf;
        else
            model.lb.lb.gain = (model.lb.lb.val - model.lb.lb.list(end))/abs(model.lb.lb.list(end));
        end
        model.lb.lb.gainlist = [model.lb.lb.gainlist model.lb.lb.gain];
        model.lb.lb.list     = [model.lb.lb.list model.lb.lb.curlist];
        model.lb.lb.it       = [model.lb.lb.it   model.lb.lb.curit];
        model.lb.lb.curlist  = [];
        model.lb.lb.curit    = [];
        if model.lb.lb.gain < 0
            strwarn = '!! dropped';
        else
            strwarn = '';
        end
        fprintf('%10s | %10s | %6.3e\n', 'Gain', strwarn, model.lb.lb.gain);
        return
    end

    % Initialise
    % ----------
    if ~isfield(model.lb, 'lb')
        model.lb.lb.val     = -inf;
        model.lb.lb.list    = []; % LB values from previous iterations
        model.lb.lb.it      = []; % x-axis values from previous iterations
        model.lb.lb.curlist = []; % LB values from current iteration
        model.lb.lb.curit   = []; % x-axis values from current iteration
    end
    
    % Compute current value
    % ---------------------
    preval = model.lb.lb.val;
    vars = fieldnames(model.lb);
    model.lb.lb.val = 0;
    for i=1:numel(vars)
        var = vars{i};
        if ~strcmpi(var, 'lb')
            model.lb.lb.val = model.lb.lb.val + model.lb.(var).val;
            if ~isfield(model.lb.(var), 'list')
                model.lb.(var).list = [];
            end
            model.lb.(var).list = [model.lb.(var).list model.lb.(var).val];
        end
    end
    model.lb.lb.curlist = [model.lb.lb.curlist model.lb.lb.val];
    diff = model.lb.lb.val - preval;

    % Update iterations (for x axis)
    % ------------------------------
    if isempty(model.lb.lb.it)
        it = 0;
    else
        it = floor(model.lb.lb.it(end)) + 1;
    end
    N  = length(model.lb.lb.curlist);
    model.lb.lb.curit = it + (0:(N-1))/N;
    
    
    % Print stuff
    % -----------
    fprintf('%10s | ', 'LB');
    if diff > 0
        fprintf('%10s | ', '(+)');
    elseif diff < 0
        fprintf('%10s | ', '(-)');
    else
        fprintf('%10s | ', '(=)');
    end
    fprintf(' %6g\n', model.lb.lb.val);
end