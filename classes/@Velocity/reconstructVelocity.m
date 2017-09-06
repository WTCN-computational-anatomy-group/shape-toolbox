function v = reconstructVelocity(obj, a1, a2, a3)
% FORMAT (v) = obj.reconstructVelocity((z), (W), (r))
% (z) - Latent coordinates [default: obj.z]
% (W) - Principal subspace [default: obj.w]
% (r) - Residual velocity field [default: obj.r]
% (v) - Output velocity [if no argout: write to obj.v]
%
% Reconstruct the velocity field from its latent PC coordinates.

    % --- Check if nothing to do
    if nargin == 1 && obj.checkarray('v')
        v = obj.v;
        return
    end
    
    % --- Default argument
    switch nargin
        case 1
            z = obj.z;
            w = obj.w;
            r = obj.r;
        case 2
            z = [];
            w = [];
            r = a1;
        case 3
            z = a1;
            w = a2;
            r = [];
        case 4
            z = a1;
            w = a2;
            r = a3;
    end
    
    % --- If all empty, get out
    if ~obj.checkarray(z) && ~obj.checkarray(r)
        if obj.Debug
            warning('Cannot reconstruct velocity: missing arrays')
        end
        v = single([]);
        return
    end
    if obj.Debug, fprintf('* reconstructVelocity\n'); end;
    
    % --- Compute linear combination
    if obj.checkarray(z)
        v = lat2vel(z, w);
    end

    % --- Add residual field
    if obj.checkarray(r)
        if obj.checkarray(z)
            v(:) = v(:) + r(:);
        else
            v = single(numeric(r));
        end
    end

    if nargout == 0
        % --- Write on disk
        obj.v.dim = size(v);
        obj.v(:)  = v(:);
        v         = obj.v;
        % --- Up-to-date
        obj.statusChanged('v');
        obj.utd.v       = true;
    end
end