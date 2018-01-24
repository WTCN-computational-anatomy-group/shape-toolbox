function setpath(component)
% FORMAT setpath((component))
% component - 'pgra'/'pgva'/['all']
%
% Set path (i.e., source dependencies) for the specified model.

    if nargin == 0
        component = {'pgva', 'pgra', 'pg', 'affine', 'diffeo'};
    end
    if ~iscell(component)
        component = {component};
    end
    
    setonepath('common');
    for i=1:numel(component)
        setonepath(component{i});
    end

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
        case 'pgva'
            addpath(fullfile(path, 'scripts'));
            addpath(fullfile(path, 'scripts', 'pgva'));
            addpath(fullfile(path, 'core', 'register'));
            addpath(fullfile(path, 'core', 'shape'));
        case 'pgra'
            addpath(fullfile(path, 'scripts'));
            addpath(fullfile(path, 'scripts', 'pgra'));
            addpath(fullfile(path, 'core', 'register'));
            addpath(fullfile(path, 'core', 'shape'));
        otherwise
            warning('Component %s not implemented yet', component)
    end
end