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
        type = 'optim';
    end
    
    ondisk = struct;
    switch lower(type)
        case 'optim'
            % Model
            ondisk.model.mu   = true;
            ondisk.model.gmu  = true;
            ondisk.model.a    = true;
            ondisk.model.w    = true;
            ondisk.model.dw   = true;
            ondisk.model.gw   = true;
            ondisk.model.hw   = true;
            ondisk.model.ww   = false;
            ondisk.model.Az   = false;
            ondisk.model.z    = false;
            ondisk.model.zz   = false;
            ondisk.model.Sz   = false;
            % Subjects
            ondisk.dat.wmu  = true;
            ondisk.dat.iphi = false;
            ondisk.dat.ipsi = false;
            ondisk.dat.v    = true;
            ondisk.dat.m    = true;
            ondisk.dat.pf   = true;
            ondisk.dat.c    = true;
            ondisk.dat.gv   = true;
            ondisk.dat.hv   = true;
            ondisk.dat.r    = true;
            ondisk.dat.gr   = false;
            ondisk.dat.hr   = false;
            ondisk.dat.z    = false;
            ondisk.dat.zz   = false;
            ondisk.dat.Sz   = false;
        case 'imageondisk'
            % Model
            ondisk.model.mu   = true;
            ondisk.model.gmu  = true;
            ondisk.model.a    = true;
            ondisk.model.w    = true;
            ondisk.model.dw   = true;
            ondisk.model.gw   = true;
            ondisk.model.hw   = true;
            ondisk.model.ww   = false;
            ondisk.model.Az   = false;
            ondisk.model.z    = false;
            ondisk.model.zz   = false;
            ondisk.model.Sz   = false;
            % Subjects
            ondisk.dat.wmu  = true;
            ondisk.dat.iphi = true;
            ondisk.dat.ipsi = true;
            ondisk.dat.v    = true;
            ondisk.dat.m    = true;
            ondisk.dat.pf   = true;
            ondisk.dat.c    = true;
            ondisk.dat.gv   = true;
            ondisk.dat.hv   = true;
            ondisk.dat.r    = true;
            ondisk.dat.gr   = true;
            ondisk.dat.hr   = true;
            ondisk.dat.z    = false;
            ondisk.dat.zz   = false;
            ondisk.dat.Sz   = false;
        case 'obswondisk'
            % Model
            ondisk.model.mu   = true;
            ondisk.model.gmu  = true;
            ondisk.model.a    = true;
            ondisk.model.w    = true;
            ondisk.model.dw   = true;
            ondisk.model.gw   = true;
            ondisk.model.hw   = true;
            ondisk.model.ww   = false;
            ondisk.model.Az   = false;
            ondisk.model.z    = false;
            ondisk.model.zz   = false;
            ondisk.model.Sz   = false;
            % Subjects
            ondisk.dat.wmu  = true;
            ondisk.dat.iphi = false;
            ondisk.dat.ipsi = false;
            ondisk.dat.v    = false;
            ondisk.dat.m    = false;
            ondisk.dat.pf   = true;
            ondisk.dat.c    = true;
            ondisk.dat.gv   = false;
            ondisk.dat.hv   = false;
            ondisk.dat.r    = false;
            ondisk.dat.gr   = false;
            ondisk.dat.hr   = false;
            ondisk.dat.z    = false;
            ondisk.dat.zz   = false;
            ondisk.dat.Sz   = false;
        case 'wondisk'
            % Model
            ondisk.model.mu   = false;
            ondisk.model.gmu  = false;
            ondisk.model.a    = false;
            ondisk.model.w    = true;
            ondisk.model.dw   = true;
            ondisk.model.gw   = true;
            ondisk.model.hw   = true;
            ondisk.model.ww   = false;
            ondisk.model.Az   = false;
            ondisk.model.z    = false;
            ondisk.model.zz   = false;
            ondisk.model.Sz   = false;
            % Subjects
            ondisk.dat.wmu  = false;
            ondisk.dat.iphi = false;
            ondisk.dat.ipsi = false;
            ondisk.dat.v    = false;
            ondisk.dat.m    = false;
            ondisk.dat.pf   = false;
            ondisk.dat.c    = false;
            ondisk.dat.gv   = false;
            ondisk.dat.hv   = false;
            ondisk.dat.r    = false;
            ondisk.dat.gr   = false;
            ondisk.dat.hr   = false;
            ondisk.dat.z    = false;
            ondisk.dat.zz   = false;
            ondisk.dat.Sz   = false;
        case 'allonram'
            % Model
            ondisk.model.mu   = false;
            ondisk.model.gmu  = false;
            ondisk.model.a    = false;
            ondisk.model.w    = false;
            ondisk.model.dw   = false;
            ondisk.model.gw   = false;
            ondisk.model.hw   = false;
            ondisk.model.ww   = false;
            ondisk.model.Az   = false;
            ondisk.model.z    = false;
            ondisk.model.zz   = false;
            ondisk.model.Sz   = false;
            % Subjects
            ondisk.dat.wmu  = false;
            ondisk.dat.iphi = false;
            ondisk.dat.ipsi = false;
            ondisk.dat.v    = false;
            ondisk.dat.m    = false;
            ondisk.dat.pf   = false;
            ondisk.dat.c    = false;
            ondisk.dat.gv   = false;
            ondisk.dat.hv   = false;
            ondisk.dat.r    = false;
            ondisk.dat.gr   = false;
            ondisk.dat.hr   = false;
            ondisk.dat.z    = false;
            ondisk.dat.zz   = false;
            ondisk.dat.Sz   = false;
        case 'allondisk'
            % Model
            ondisk.model.mu   = true;
            ondisk.model.gmu  = true;
            ondisk.model.a    = true;
            ondisk.model.w    = true;
            ondisk.model.dw   = true;
            ondisk.model.gw   = true;
            ondisk.model.hw   = true;
            ondisk.model.ww   = true;
            ondisk.model.Az   = true;
            ondisk.model.z    = true;
            ondisk.model.zz   = true;
            ondisk.model.Sz   = true;
            % Subjects
            ondisk.dat.wmu  = true;
            ondisk.dat.iphi = true;
            ondisk.dat.ipsi = true;
            ondisk.dat.v    = true;
            ondisk.dat.m    = true;
            ondisk.dat.pf   = true;
            ondisk.dat.c    = true;
            ondisk.dat.gv   = true;
            ondisk.dat.hv   = true;
            ondisk.dat.r    = true;
            ondisk.dat.gr   = true;
            ondisk.dat.hr   = true;
            ondisk.dat.z    = true;
            ondisk.dat.zz   = true;
            ondisk.dat.Sz   = true;
    end
    
end