function swz = uncertaintyWZ(obj)

    if nargin == 1 && obj.checkarray('swz')
        swz = obj.swz;
        return
    end
    
    obj.uncertaintyZ();
    if ~obj.checkarray('sz') || ~obj.checkarray('w')
        if obj.Debug
            warning('Cannot compute WZ''s uncertainty. Missing arrays.')
        end
        return
    end
    if obj.Debug, fprintf('* uncertaintyWZ\n'); end;
    
    sz = obj.sz;
    w  = obj.w;
    
    swz = latentToTensorCov(w, sz);
    
    obj.swz.dim = size(swz);
    obj.swz(:)  = swz(:);
    obj.statusChanged('swz');
    swz         = obj.swz;
    
end