function out = mapProperties(osource, psource, otarget, ptarget, ind, f)
% FORMAT mapProperties(osource, psource, otarget, ptarget, ind, f)
% osource - If object:      handle source object
%           If 'value':     psource must be a value
%           If 'values':    psource must be a list of values
% psource - If char:        source's property to map.
%           Else:           value(s) to map
% otarget - If char: property containing the list of target handle objects.
%           If cell: list of target handle objects.
% ptarget - Target's property to set [default: psource].
% ind     - Indices of the targets to set [default: all].
% f       - Function to apply to the mapped values [default: identity].
%
% Map a property value to a series of objects.
% > mapProperties(obj, 'MatchingTerm', 'Velocities', [1 3 5]);
% > mapProperties('value', 1, {v1, v2, v3}, 'MaxGNIt', 1:3, @(X) 2*X);
% > mapProperties('values', {1, 2, 3}, {v1, v2, v3}, 'MaxGNIt');

    % --- Default arguments
    if nargin < 6
        f = @(X) X;
        if nargin < 5
            ind = 1:length(osource.(otarget));
            if nargin < 4
                ptarget = psource;
            end
        end
    end
    
    % --- Prepare value to map
    one_value = true;
    if ichar(osource)
        if strcmpi(osource, 'value')
            value = f(psource);
        elseif strcmpi(osource, 'values')
            value = cellfun(f, psource, 'UniformOutput', false);
            one_value = false;
        end
    else
        st = osource.disableListeners(ptarget);
        try
            value = f(osource.(psource));
        catch e
            warning('mapProperty failed');
            osource.enableListeners(st, ptarget);
            rethrow(e)
        end
        osource.enableListeners(st, ptarget);
    end
    
    % --- Map
    if ischar(otarget)
        otarget = osource.(otarget);
    end
    if one_value
        for i=ind
            otarget{i}(ptarget) = value;
        end
    else
        for i=ind
            otarget{i}(ptarget) = value{i};
        end
    end
    
    % --- Output
    if ichar(osource)
        out = otarget;
    else
        out = osource;
    end
    
end