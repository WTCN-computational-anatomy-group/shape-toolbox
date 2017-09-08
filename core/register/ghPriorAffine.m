function [g,h] = ghPriorAffine(q, reg, varargin)
% FORMAT [g,(h)] = obj.ghPriorAffine(q, reg, ...)
%
% ** Required **
% q   - Current affine parameters
% reg - Precision (or 'inverse covariance' or 'regularisation') matrix.
% ** Output **
% g  - Gradient ([nq 1])
% h  - Hessian  ([nq nq])
%
% - Compute gradient/hessian of the Affine prior term w.r.t. parameters.
% - The affine basis is supposed to be ordered so that the non regularized 
%   parameters are first (and the regularised ones last). Usually,
%   translations and rotations are not regularised.

    % --- Parse inputs
    p = inputParser;
    p.FunctionName = 'ghPriorAffine';
    p.addRequired('q',     @checkarray);
    p.addRequired('reg',   @checkarray);
    p.addParameter('output',   []);
    p.addParameter('debug',    false,   @isscalar);
    p.parse(q, reg, varargin{:});
    output  = p.Results.output;
    debug   = p.Results.debug;
    
    if debug, fprintf('* ghPriorAffine\n'); end;
    
    
    % --- Dim info
    nq = size(q, 1);
    npq = size(reg, 1);
    
    % --- Load data
    q    = numeric(q);
    reg  = numeric(reg);
    % The first basis (rigid) have no prior
    sub  = (nq-npq+1):nq; 
    
    % --- Compute grad and hess
    g = zeros([nq 1]);
    g(sub) = reg * q(sub);
    if nargout > 1
        h = zeros(nq);
        h(sub,sub) = reg;
    end
    
    % --- Write on disk
    if ~iscell(output)
        output = {output};
    end
    if numel(output) < 2
        output = [output {[]}];
    end
    if ~isempty(output{1})
        g = saveOnDisk(g, output{1}, 'name', 'g');
    end
    if nargout > 1 && ~isempty(output{2})
        h = saveOnDisk(h, output{2}, 'name', 'h');
    end
    
end