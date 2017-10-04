function ondisk = default_ondisk(type)
% FORMAT ondisk = default_ondisk(type)
% type - On/Off disk scheme:
%        * 'ImageOnDisk': All "images" (i.e., defined on a lattice) are on
%                         disk, the rest is on RAM.
%        * 'ObsWOnDisk':  Subspace and "observed images" (f, mu, pf...) 
%                         are on disk, the rest is on RAM.
%        * 'WOnDisk':     Only the subspace (w, gw, hw, dw) is on disk.
%        * 'AllOnRAM':    Everything on RAM.
%        * 'AllOnDisk':   Everything on disk.

    if nargin < 1
        type = 'imageondisk';
    end
    
    ondisk = struct;
    switch lower(type)
        case 'imageondisk'
            % Model
            ondisk.model.mu   = true;
            ondisk.model.gmu  = true;
            ondisk.model.a    = true;
            ondisk.model.w    = true;
            ondisk.model.dw   = true;
            ondisk.model.g    = true;
            ondisk.model.h    = true;
            ondisk.model.ww   = false;
            ondisk.model.A    = false;
            ondisk.model.A0   = false;
            ondisk.model.z    = false;
            ondisk.model.zz   = false;
            ondisk.model.S    = false;
            % Subjects
            ondisk.dat.wmu  = true;
            ondisk.dat.iphi = true;
            ondisk.dat.ipsi = true;
            ondisk.dat.v    = true;
            ondisk.dat.pf   = true;
            ondisk.dat.c    = true;
            ondisk.dat.gv   = true;
            ondisk.dat.hv   = true;
            ondisk.dat.z    = false;
            ondisk.dat.zz   = false;
            ondisk.dat.S    = false;
        case 'obswondisk'
            % Model
            ondisk.model.mu   = true;
            ondisk.model.gmu  = true;
            ondisk.model.a    = true;
            ondisk.model.w    = true;
            ondisk.model.dw   = true;
            ondisk.model.g    = true;
            ondisk.model.h    = true;
            ondisk.model.ww   = false;
            ondisk.model.A    = false;
            ondisk.model.A0   = false;
            ondisk.model.z    = false;
            ondisk.model.zz   = false;
            ondisk.model.S    = false;
            % Subjects
            ondisk.dat.wmu  = true;
            ondisk.dat.iphi = false;
            ondisk.dat.ipsi = false;
            ondisk.dat.v    = false;
            ondisk.dat.pf   = true;
            ondisk.dat.c    = true;
            ondisk.dat.gv   = false;
            ondisk.dat.hv   = false;
            ondisk.dat.z    = false;
            ondisk.dat.zz   = false;
            ondisk.dat.S    = false;
        case 'wondisk'
            % Model
            ondisk.model.mu   = false;
            ondisk.model.gmu  = false;
            ondisk.model.a    = false;
            ondisk.model.w    = true;
            ondisk.model.dw   = true;
            ondisk.model.g    = true;
            ondisk.model.h    = true;
            ondisk.model.ww   = false;
            ondisk.model.A    = false;
            ondisk.model.A0   = false;
            ondisk.model.z    = false;
            ondisk.model.zz   = false;
            ondisk.model.S    = false;
            % Subjects
            ondisk.dat.wmu  = false;
            ondisk.dat.iphi = false;
            ondisk.dat.ipsi = false;
            ondisk.dat.v    = false;
            ondisk.dat.pf   = false;
            ondisk.dat.c    = false;
            ondisk.dat.gv   = false;
            ondisk.dat.hv   = false;
            ondisk.dat.z    = false;
            ondisk.dat.zz   = false;
            ondisk.dat.S    = false;
        case 'allonram'
            % Model
            ondisk.model.mu   = false;
            ondisk.model.gmu  = false;
            ondisk.model.a    = false;
            ondisk.model.w    = false;
            ondisk.model.dw   = false;
            ondisk.model.g    = false;
            ondisk.model.h    = false;
            ondisk.model.ww   = false;
            ondisk.model.A    = false;
            ondisk.model.A0   = false;
            ondisk.model.z    = false;
            ondisk.model.zz   = false;
            ondisk.model.S    = false;
            % Subjects
            ondisk.dat.wmu  = false;
            ondisk.dat.iphi = false;
            ondisk.dat.ipsi = false;
            ondisk.dat.v    = false;
            ondisk.dat.pf   = false;
            ondisk.dat.c    = false;
            ondisk.dat.gv   = false;
            ondisk.dat.hv   = false;
            ondisk.dat.z    = false;
            ondisk.dat.zz   = false;
            ondisk.dat.S    = false;
        case 'allondisk'
            % Model
            ondisk.model.mu   = true;
            ondisk.model.gmu  = true;
            ondisk.model.a    = true;
            ondisk.model.w    = true;
            ondisk.model.dw   = true;
            ondisk.model.g    = true;
            ondisk.model.h    = true;
            ondisk.model.ww   = true;
            ondisk.model.A    = true;
            ondisk.model.A0   = true;
            ondisk.model.z    = true;
            ondisk.model.zz   = true;
            ondisk.model.S    = true;
            % Subjects
            ondisk.dat.wmu  = true;
            ondisk.dat.iphi = true;
            ondisk.dat.ipsi = true;
            ondisk.dat.v    = true;
            ondisk.dat.pf   = true;
            ondisk.dat.c    = true;
            ondisk.dat.gv   = true;
            ondisk.dat.hv   = true;
            ondisk.dat.z    = true;
            ondisk.dat.zz   = true;
            ondisk.dat.S    = true;
    end
    
end