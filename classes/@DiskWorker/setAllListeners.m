function setAllListeners(obj)
% FORMAT obj.setAllListeners()
%
% Creates and stores "classical" listeners on public niftis.
% Three listeners are set for each:
% * PreSet:  Saves the previous value in a private structure.
% * PostSet: Uses the input and the old value to decide how to set the new
%            value.
% * PreGet:  Call a dedicated updater to insure that the returned array is
%            adequately set.
% All these listeners are stored in obj.alllisteners (a complete list)
% and obj.listeners.[public_name] (a fast access map).

    function change_status(obj, name)
        obj.statusChanged(name);
        obj.utd.(name) = true;
    end

    obj.alllisteners = {};
    obj.listeners = struct;
    
    % --- GetObservable properties
    p = findAttrValue(obj, 'GetObservable', 'true');
    for i=1:numel(p)
        % - set up-to-date if needed
        if isa(obj.(p{i}), 'nifti') && isfield(obj.nii2dat, p{i})
            utd_name = obj.nii2dat.(p{i});
        else
            utd_name = p{i};
        end
        if ~isfield(obj.utd, utd_name)
            obj.utd.(utd_name) = false;
        end
        % - init listener struct if needed
        if ~isfield(obj.listeners, p{i})
            obj.listeners.(p{i}) = {};
        end
        % - Register listener
        obj.listeners.(p{i}){end+1} = addlistener(obj, p{i}, 'PreGet',  @obj.getProperty);  % Check up-to-date
        obj.alllisteners{end+1} = obj.listeners.(p{i}){end};
        obj.alllisteners{end}.Enabled = false;
    end
    
    % --- SetObservable properties
    p = findAttrValue(obj, 'SetObservable', 'true');
    for i=1:numel(p)
        % - init listener struct if needed
        if ~isfield(obj.listeners, p{i})
            obj.listeners.(p{i}) = {};
        end
        % - Register listener
        if isa(obj.(p{i}), 'nifti')
            obj.listeners.(p{i}){end+1} = addlistener(obj, p{i}, 'PreSet',  @obj.saveArray);    % Save previous data
            obj.alllisteners{end+1} = obj.listeners.(p{i}){end};
            obj.alllisteners{end}.Enabled = false;
            obj.listeners.(p{i}){end+1} = addlistener(obj, p{i}, 'PostSet', @obj.setArray);     % Call appropriate setter
            obj.alllisteners{end+1} = obj.listeners.(p{i}){end};
            obj.alllisteners{end}.Enabled = false;
        else
            obj.listeners.(p{i}) = {addlistener(obj, p{i}, 'PostSet',  @(s,e) change_status(obj, s.Name) )};  % Change up-to-date status
            obj.alllisteners{end+1} = obj.listeners.(p{i}){1};
            obj.alllisteners{end}.Enabled = false;
        end
    end
    
    % --- Enable everything
    obj.enableListeners();
end