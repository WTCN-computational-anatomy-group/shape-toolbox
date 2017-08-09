function swz = uncertaintyWZ(obj)

    if ~obj.checkarray('sz')
        obj.uncertaintyZ();
    end
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
    obj.utd.swz = true;
    swz         = obj.swz;
    
end