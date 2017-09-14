function [g,h] = ghPriorLatent(z, reg, varargin)
% FORMAT [g,(h)] = obj.ghPriorLatent(z, reg, ...)
%
% ** Required **
% z   - Current latent parameters
% reg - Precision (or 'inverse covariance' or 'regularisation') matrix.
% ** Output **
% g  - Gradient ([nk 1])
% h  - Hessian  ([nk nk])
%
% - Compute gradient/hessian of the latent prior term w.r.t. parameters.

    % --- Parse inputs
    p = inputParser;
    p.FunctionName = 'ghPriorLatent';
    p.addRequired('z',     @checkarray);
    p.addRequired('reg',   @checkarray);
    p.addParameter('output',   []);
    p.addParameter('debug',    false,   @isscalar);
    p.parse(z, reg, varargin{:});
    output  = p.Results.output;
    debug   = p.Results.debug;
    
    if debug, fprintf('* ghPriorLatent\n'); end;
    
    % --- Load data
    z    = numeric(z);
    reg  = numeric(reg);
    
    % --- Compute grad and hess
    g = reg * z;
    if nargout > 1
        h = reg;
    end

    % --- Write on disk
    if ~iscell(output)
        output = {output};
    end
    if numel(output) < 2
        output = [output {[]}];
    end
    if ~isempty(output{1})
        g = saveOnDisk(output{1}, g, 'name', 'g');
    end
    if nargout > 1 && ~isempty(output{2})
        h = saveOnDisk(output{2}, h, 'name', 'h');
    end
end