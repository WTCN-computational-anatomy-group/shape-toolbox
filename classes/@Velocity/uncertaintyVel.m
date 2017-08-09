function sv = uncertaintyVel(obj)

    if ~obj.checkarray('sr')
        obj.uncertaintyR();
    end
    if ~obj.checkarray('swz')
        obj.uncertaintyWZ();
    end

    if ~obj.checkarray('swz') && ~obj.checkarray('sr')
        if obj.Debug
            warning('Cannot compute velocity''s uncertainty. Missing arrays.')
        end
        return
    end
    if obj.Debug, fprintf('* uncertaintyVel\n'); end;
    
    if obj.checkarray('swz')
        swz = numeric(obj.swz);
    else
        swz = [];
    end
    
    if obj.checkarray('sr')
        sr = numeric(obj.sr);
    else
        sr = [];
    end
    
    if ~isempty(swz)
        if ~isempty(sr)
            sv = swz + obj.SigmaR^2 * sr;
        else
            sv = swz;
        end
    elseif ~isempty(sr)
        sv = obj.SigmaR^2 * sr;
    else
        sv = [];
    end
    
    obj.sv.dim = size(sv);
    obj.sv(:) = sv(:);
    obj.statusChanged('sv');
    obj.utd.sv = true;
    sv = obj.sv;
    
end