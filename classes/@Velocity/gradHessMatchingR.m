% //!\\ Not used at the moment.
% The current implementation call gradHessMatchingVel directly from
% computeGradHessR.

function [g, h] = gradHessMatchingR(obj, g, h)
% FORMAT [g, h] = obj.gradHessMatchingZ((g), (h))
% obj - A velocity object
% (g) - Gradient of the objective function w.r.t. residual field.
%       Incremented by this function.
% (h) - Hessian of the objective function w.r.t. residual field.
%       Incremented by this function.
%
% Object fields [mu, pf, pvox, MatchingTerm] are used by this function 
% and should thus be correctly set before call.
%
% Add the matching term part to Gradient and Hessian
% This function should only be called by computeGradHess.

    if nargin < 3
        h = [];
        if nargin < 2
            g = [];
        end
    end

    % --- Check all arrays are ready to be used
    obj.reconstructTemplate();
    obj.computeTemplateGrad();
    obj.pushImage();
    if ~obj.checkarray('w') || ~obj.checkarray('mu') || ...
       ~obj.checkarray('pf')  || ~obj.checkarray('pvox') || ...
       ~obj.checkarray('gmu')
        warning('Cannot compute gradient and hessian: missing arrays\n');
        return
    end
    if obj.Debug, fprintf('* gradHessMatchingR\n');  end;

    % --- Read data
    dim         = [size(obj.w) 1 1 1];
    dim         = dim(1:5);
    dim_lattice = dim(1:3);
    dim_vector  = dim(4);
    
    if nargin < 3
        h = zeros([dim_lattice dim_vector*(dim_vector+1)/2]);
        if nargin < 2
            g = zeros([dim_lattice dim_vector]);
        end
    end

    % Loop on Z-slices to save memory
    for x3=1:dim_lattice(3)

        % --- Load template and pushed image
        gmu = single(obj.gmu(:,:,x3,:,:));   % Template spatial gradients
        mu  = single(obj.mu(:,:,x3,:));      % Template
        f   = single(obj.pf(:,:,x3,:));      % Pushed image
        c   = single(obj.pvox(:,:,x3,:));    % Pushed voxel count

        % --- Compute grad/hess w.r.t. the complete velocity
        [g1, h1] = gradHessMatchingVel(obj, mu, f, c, gmu);
        clear mu f c gmu
        
        % --- Increment grad/hess
        g = g + g1;
        h = h + h1;
        clear g1 h1

    end % End loop on Z-slices
    
    % --- Multiply by the ammount of noise (that way ?)
    g = g * obj.SigmaR;
    h = h * (obj.SigmaR^2);
    
    % Just in case
    g(~isfinite(g)) = 0;
    h(~isfinite(h)) = 0;
end