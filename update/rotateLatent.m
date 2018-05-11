function dat = rotateLatent(R, dat, ~, ~)
% FORMAT dat = rotateLatent(dat, model, opt)
% dat   - Subject-specific data
% model - Model-specific data
% opt   - Options
%
% All inputs are structures that can either be in memory or on disk in the
% form of a mat file. In the latter case, it is read and, if needed,
% written back.
%--------------------------------------------------------------------------
% Rotate latent shape coordinates
%__________________________________________________________________________
% Copyright (C) 2018 Wellcome Centre for Human Neuroimaging

    % =====================================================================
    % Read input from disk (if needed)
    [dat, datpath] = fileToStruct(dat);
   
    % =====================================================================
    % Rotate
    z        = R * dat.z.z(:);
    dat.z.z  = z;
    dat.z.zz = z*z';
    dat.z.S  = R * dat.z.S * R';
                             
    % =====================================================================
    % Exit
    dat = structToFile(dat, datpath);
end