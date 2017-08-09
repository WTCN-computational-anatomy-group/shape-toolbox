classdef DiskWorker < handle
% CLASS DiskWorker
%
% Base class for processing classes that intend to memory map data.
% It defines utilities for setting niftis, listening to changes, keeping
% track of dependencies and array states...
%
% For each nifti property (we advise to use CaptializedNames), a dependent
% file_array property (we advise to use lowercase_names) should be define, 
% with the folloowing getter and setter:
%   function val = get.private_name(obj),
%       s = obj.disableListeners('PublicName');
%       val = obj.PublicName.dat;
%       obj.enableListeners(s, 'PublicName');
%   end
%   function set.private_name(obj, val)
%       s = obj.disableListeners('PublicName');
%       obj.PublicName.dat = val;
%       obj.enableListeners(s, 'PublicName');
%   end
% You can also use
% > value = obj.getDat('private_name')
% and
% > obj.setDat('private_name', value)
% which are equivalent and already defined.
%
% Child classes should define
% - nii2dat
%       A constant structure defining a translation between public nifti
%       (i.e., properties which value is a nifti object) and private file
%       arrays (i.e., dependent properties which redirect towards the 
%       nifti's file_array).
% - dependencies
%       A constant structure defining dependencies between arrays or 
%       properties. It should be of the form:
%           s.d1 = {'u1', 'u2', 'u3'},
%       meaning that computing the array or property d1 requires arrays or 
%       properties u1, u2 and u3. For arrays, private names should be used.
% - updaters
%       A constant structure associating an update function with each
%       (private or public) property: It should be of the form:
%           s.d = @computeD
%
% Finally, at the end of the constructor, call obj.setAllListeners().


    %% == Abstract properties =============================================
    properties (Abstract, Constant, Access = protected)
        dependencies        % Structure describing upstream dependencies between "stuff": D depends on [U1, ..., Un]
        nii2dat             % Translation between public property names and the ones used in utd/dependencies
        updaters            % Functions to call to update a specific property
        % --- Auto properties
        dat2nii             % Just call invertTranslation(nii2dat) to initialise it
        downdep             % Just call invertDependencies(dependencies) to initialise it
    end
    %% == Listeners utils =================================================
    properties (Access = protected)
        % --- Up-To-Date info ---
        % We defined PreGet listeners on arrays to check if they are
        % up-to-date before returning them. These variables allow us to
        % know if they are.
        
        utd           = struct()    % Structure where fields are array public names
        private                     % Temporary store data to help listeners
        listeners     = struct()    % References to listeners
        alllisteners  = {}          % References to listeners
    end
    methods (Access = protected)
        setAllListeners(obj)
        rmAllListeners(obj)
        prev_state = disableListeners(obj, varargin)
        enableListeners(obj, varargin)
        saveArray(obj, s, ~)
        setArray(obj, s, ~)
        getProperty(obj, s, ~)
        statusChanged(obj, varargin)
        ok = checkarray(obj, array)
    end
    methods (Static, Access = protected)
       
        setDat(obj, prop, value)
        value = getDat(obj, prop)
        
        function s = invertTranslation(nii2dat)
            s = struct;
            values = fieldnames(nii2dat);
            for i=1:numel(values)
                s.(nii2dat.(values{i})) = values{i};
            end
        end
        
        
        function s = invertDependencies(dep)
            s = struct;
            fields = fieldnames(dep);
            for i=1:numel(fields)
                parents = dep.(fields{i});
                for j=1:numel(parents)
                    if ~isfield(s, parents{j})
                        s.(parents{j}) = {};
                    end
                    s.(parents{j}){end+1} = fields{i};
                end
            end
            fields = fieldnames(s);
            for i=1:numel(fields)
                s.(fields{i}) = unique(s.(fields{i}));
            end
        end
    end
end