function [par, loop] = autoParLoop(par, loop, isfa, nz, nc)
% FORMAT [par, loop] = autoParLoop(par, loop, isfa, nz, (nc))
% par  - User parameter for the parallelisation scheme.
% loop - User parameter for the splitting scheme.
% isfa - Is the input a file_array ?
% nz   - Number of z slices
% nc   - Number of component [default: do not loop over components]
%
% - Tries do determine the optimal split and parallelisation scheme based 
%   on the context (are we already a parallelised job) and the input type
%   (is it a file_array, how many slices and components).
% - User arguments are usually kept unless they are clearly inefficient
%   (parallelising only 1 slice for example).
% - The aims are (roughly):
%   1) Minimise memory usage (split loading if the input is a file_array)
%   2) Maximise parallelisation usage
%   3) Avoid unnecessary parallelisation (over a loop of length 1)

    if nargin < 5
        nc = 1;
    end

    
    if strcmpi(loop, 'auto')
        loop = '';
    end
    
    % --- Splitting scheme
    if isempty(loop)
        if nz > 1
            loop = 'slice';
        elseif nc > 1
            loop = 'component';
        else
            loop = 'none';
        end
    end
    
    % --- Parallelisation scheme
    try
        parid = getCurrentTask();
    catch
        parid = '';
    end
    if ~isempty(parid)
        par = 0;
    else
        if islogical(par)
            if par
                par = inf;
            else
                par = 0;
            end
        end
    end
    
    % --- Update splitting scheme based on parallelisation scheme
    if par == 0 && ~isfa
        loop = 'none';
    end

end