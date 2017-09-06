function prev_state = disableListeners(obj, varargin)
% FORMAT prev_state = obj.disableListeners(('PublicName1', ...))
%
% Temporary disable some listeners based on the public name of the property
% they are linked to. If no name is provided, all listeners are disabled.
% An encoding of the listeners state before call is returned, so that it
% can be restored later.

% 2/8/17 I am removing all the cellfun and use for-loop instead as I
% discovered those are (surprisingly) faster.

    if numel(varargin) == 0
    % --- Disable all
        try
            prev_state = zeros([1 numel(obj.alllisteners)], 'logical');
        catch
            prev_state = zeros([1 numel(obj.alllisteners)], 'uint8');
        end
        for j=1:numel(obj.alllisteners)
            prev_state(j) = obj.alllisteners{j}.Enabled;
            obj.alllisteners{j}.Enabled = false;
        end
        
    else
    % --- Disable arguments
        try
            prev_state = zeros([1 numel(varargin)], 'logical');
        catch
            prev_state = zeros([1 numel(varargin)], 'uint8');
        end
        for i=1:numel(varargin)
            if isfield(obj.listeners, varargin{i})
                prev_state(i) = false;
                for j=1:numel(obj.listeners.(varargin{i}))
                    prev_state(i) = prev_state(i) || obj.listeners.(varargin{i}){j}.Enabled;
                    obj.listeners.(varargin{i}){j}.Enabled = false;
                end
            else
                prev_state(i) = false;
            end
        end
    end
end

%--------------------------------------------------------------------------
% OLD cellfun implem
%--------------------------------------------------------------------------
% function prev_state = disable(X)
%     prev_state = X.Enabled;
%     X.Enabled  = false;
% end
% if numel(varargin) == 0
%     % Disable all
%     prev_state = cellfun(@disable, obj.alllisteners);
% else
%     % Disable arguments
%     prev_state = zeros([1 numel(varargin)], 'logical');
%     for i=1:numel(varargin)
%         if isfield(obj.listeners, varargin{i})
%             prev_state(i) = any(cellfun(@disable, obj.listeners.(varargin{i})));
%         else
%             prev_state(i) = false;
%         end
%     end
% end