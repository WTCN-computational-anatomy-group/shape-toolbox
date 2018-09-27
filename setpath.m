function setpath(component)
% FORMAT setpath((component))
% component - 'pgra'/'pgva'/['all']
%
% Set path (i.e., source dependencies) for the specified model.

    warning('off','MATLAB:dispatcher:nameConflict')

    if nargin == 0
        component = {'shape_model'};
    end
    if ~iscell(component)
        component = {component};
    end
    
    setonepath('common');
    for i=1:numel(component)
        setonepath(component{i});
    end

    warning('on','MATLAB:dispatcher:nameConflict')
end

function setonepath(component)
    
    path = fileparts(which('setpath'));
    switch component
        case 'common'
            addpath(fullfile(path, 'backward-compatibility'));
            addpath(fullfile(path, 'utility-functions'));
            if exist(fullfile(path, '..', 'auxiliary-functions'), 'dir')
                addpath(fullfile(path, '..', 'auxiliary-functions'));
            end
            if exist(fullfile(path, '..', 'distributed-computing'), 'dir')
                addpath(fullfile(path, '..', 'distributed-computing'));
            end
        case 'shape_model'
            addpath(fullfile(path, 'scripts'));
            addpath(fullfile(path, 'scripts', 'shape_model'));
            addpath(fullfile(path, 'core', 'register'));
            addpath(fullfile(path, 'core', 'shape'));
            addpath(fullfile(path, 'update'));
        otherwise
            warning('Component %s not implemented yet', component)
    end
end