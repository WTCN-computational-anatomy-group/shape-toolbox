function psi = reconstructPsi(obj, A, phi, lat)

    % --- Return if nothing to do
    if nargin == 1 && checkarray('psi')
        psi = obj.ipsi;
        return
    end

    % --- Default arguments
    if nargin < 4
        lat = [size(obj.mu) 1];
        lat = lat(1:3);
        if nargin < 3
            obj.exponentiateVelocity('phi');
            phi = obj.phi;
            if nargin < 2
                obj.exponentiateAffine();
                A = obj.A;
            end
        end
    end
    
    % --- Check the everything is alright
    if ~checkarray(phi) && (~checkarray(lat) ||  ~checkarray(A))
        if obj.Debug
            warning('Cannot reconstruct psi: missing arrays');
        end
        return
    end
    if obj.Debug, fprintf('* reconstructPsi\n'); end;
    
    % --- Load data
    phi = numeric(phi);
    
    % --- Reconstruct
    if isempty(phi)
        id = warps('identity', lat);
        psi = warps('compose', obj.Mmu \ A \ obj.Mf, id);
    elseif isempty(A)
        psi = phi;
    else
        psi = warps('compose', obj.Mmu \ A \ obj.Mf, phi);
    end
    
    % --- Write
    if nargout == 0
        obj.psi.dim = size(psi);
        obj.psi(:) = psi(:);
        obj.statusChanged('psi');
    end
    
end