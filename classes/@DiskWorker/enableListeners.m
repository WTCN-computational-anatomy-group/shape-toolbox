function enableListeners(obj, varargin)
% FORMAT objenableListeners((prev_state), ('PublicName1', ...))
%
% (Re)-enable some listeners based on the public name of the property
% they are linked to. If no name is provided, all listeners are enabled.
% An encoding of the listeners state can be provided, so that it
% is restored (rather than forcing an "enabled" state).

% 2/8/17 I am removing all the cellfun and use for-loop instead as I
% discovered those are (surprisingly) faster.

    % --- Is a state provided ?
    prev_state = [];
    if numel(varargin) > 0 && islogical(varargin{1})
        prev_state = varargin{1};
        varargin = varargin(2:end);
    end
    
    if numel(varargin) == 0
    % --- Enable all
        if isempty(prev_state)
            for j=1:numel(obj.alllisteners)
                obj.alllisteners{j}.Enabled = true;
            end
        else
            for j=1:numel(obj.alllisteners)
                obj.alllisteners{j}.Enabled = prev_state(j);
            end
        end
        
    else
    % --- Enable arguments
        if isempty(prev_state)
            for i=1:numel(varargin)
                if isfield(obj.listeners, varargin{i})
                    for j=1:numel(obj.listeners.(varargin{i}))
                        obj.listeners.(varargin{i}){j}.Enabled = true;
                    end
                end
            end
        else
            for i=1:numel(varargin)
                if isfield(obj.listeners, varargin{i})
                    for j=1:numel(obj.listeners.(varargin{i}))
                        obj.listeners.(varargin{i}){j}.Enabled = prev_state(i);
                    end
                end
            end
        end
    end

end

%--------------------------------------------------------------------------
% OLD cellfun implem
%--------------------------------------------------------------------------
% function enable(X, V)
%     if nargin < 2
%         V = true;
%     end
%     X.Enabled = V;
% end
% prev_state = [];
% if numel(varargin) > 0 && islogical(varargin{1})
%     prev_state = varargin{1};
%     varargin = varargin(2:end);
% end
% if numel(varargin) == 0
%     % Enable all
%     if isempty(prev_state)
%         cellfun(@enable, obj.alllisteners);
%     else
%         cellfun(@enable, obj.alllisteners, num2cell(prev_state));
%     end
% else
%     % Enable arguments
%     if isempty(prev_state)
%         for i=1:numel(varargin)
%             if isfield(obj.listeners, varargin{i})
%                 cellfun(@enable, obj.listeners.(varargin{i}));
%             end
%         end
%     else
%         for i=1:numel(varargin)
%             if isfield(obj.listeners, varargin{i})
%                 cellfun(@enable, obj.listeners.(varargin{i}), ...
%                         {prev_state(i), prev_state(i), prev_state(i)});
%             end
%         end
%     end
% end