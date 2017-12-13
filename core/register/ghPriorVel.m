function g = ghPriorVel(v, varargin)
% FORMAT g = ghPriorVel(v, (vs), (prm), (bnd))
%
% ** Required **
% v   - Initial velocity
% ** Optional **
% vs  - Voxel size of the template lattice [1 1 1]
% prm - Parameters of the L differential operator
%       [0.0001 0.001 0.2 0.05 0.2]
% bnd - L differential operator boundary conditions (0/1/2/3) [0]
% ** Output **
% g   - Gradient
%
% Gradient of the *negative* log-likelihood of the prior term for the
% initial velocity (Lv).

    % --- Parse inputs
    p = inputParser;
    p.FunctionName = 'ghPriorVel';
    p.addRequired('v',         @checkarray);
    p.addOptional('vs',        [1 1 1], @(X) length(X) == 3);
    p.addOptional('prm',       [0.0001 0.001 0.2 0.05 0.2], @(X) length(X) == 5);
    p.addOptional('bnd',       0, @(X) isscalar(X) && isnumeric(X));
    p.addParameter('output',   []);
    p.addParameter('debug',    false,   @isscalar);
    p.parse(v, varargin{:});
    vs      = p.Results.vs;
    prm     = p.Results.prm;
    bnd     = p.Results.bnd;
    output  = p.Results.output;
    debug   = p.Results.debug;
    
    if debug, fprintf('* ghPriorVel\n'); end
    
    spm_diffeo('boundary', bnd);
    g = spm_diffeo('vel2mom', single(numeric(v)), double([vs prm]));
    
    if ~isempty(output)
        g = writeOnDisk(output, g, 'name', 'g');
    end

end