function ipsi = reconstructIPsi(obj, A, iphi, lat)

    % --- Return if nothing to do
    if nargin == 1 && obj.checkarray('ipsi')
        ipsi = obj.ipsi;
        return
    end

    % --- Default arguments
    if nargin < 4
        lat = [size(obj.f) 1];
        lat = lat(1:3);
        if nargin < 3
            obj.exponentiateVelocity('iphi');
            iphi = obj.iphi;
            if nargin < 2
                obj.exponentiateAffine();
                A = obj.A;
            end
        end
    end
    
    % --- Check the everything is alright
    if  ~obj.checkarray(lat) || (~obj.checkarray(A) && ~obj.checkarray(iphi))
        if obj.Debug
            warning('Cannot reconstruct ipsi: missing arrays');
        end
        return
    end
    if obj.Debug, fprintf('* reconstructIPsi\n'); end;
    
    % --- Load data
    iphi = numeric(iphi);
    
    % --- Reconstruct
    id = warps('identity', lat);
    if isempty(iphi)
        ipsi = warps('compose', obj.Mmu \ A \ obj.Mf, id);
    elseif isempty(A)
        ipsi = warps('compose', iphi, id);
    else
        ipsi = warps('compose', iphi, obj.Mmu \ A \ obj.Mf, id);
    end
    
    % --- Write
    if nargout == 0
        obj.ipsi.dim = size(ipsi);
        obj.ipsi(:) = ipsi(:);
        obj.statusChanged('ipsi');
    end
    
end