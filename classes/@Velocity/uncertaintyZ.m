function sz = uncertaintyZ(obj)

    if nargin == 1 && obj.checkarray('sz')
        sz = obj.sz;
        return
    end

    if ~obj.checkarray('hz')
        if obj.Debug
            warning('Cannot compute Z''s uncertainty. Missing arrays.')
        end
        return
    end
    if obj.Debug, fprintf('* uncertaintyZ\n'); end;
    
    hz = obj.hz;
    
    obj.sz.dim = size(hz);
    sz = inv(numeric(hz));
    obj.sz(:) = sz(:);
    obj.statusChanged('sz');
    sz = obj.sz;
    
end