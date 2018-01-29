function [pf, c, bb] = pushImage(ipsi, f, varargin)
% _________________________________________________________________________
%
% Push the image to template space
%
% -------------------------------------------------------------------------
%
% FORMAT ([pf, (c), (bb)]) = pushImage(ipsi, f, (lat), ...)
%
% REQUIRED
% --------
% ipsi - Inverse transform (warps mu to f).
% f    - Observed image.
%
% OPTIONAL
% --------
% lat  - Output lattice [same as input]
%
% KEYWORD ARGUMENTS
% -----------------
% loop   - How to split processing ('none', 'component') [auto]
% par    - If true, parallelise processing over classes/modalities [false]
% output - Cell of file_array to use as output [return numeric array]
% debug  - Print debugging info [false]
%
% OUTPUT
% ------
% pf   - Pushed image in template space 
% c    - Pushed voxel count
% bb   - Bounding box: correspondance between lat indices and the saved
%        volume
%
% -------------------------------------------------------------------------
% - Interpolation is trilinear.
% - Boundary conditions are circulant.
% _________________________________________________________________________

    % --- Parse inputs
    p = inputParser;
    p.FunctionName = 'pushImage';
    p.addRequired('ipsi', @checkarray);
    p.addRequired('f',    @checkarray);
    p.addOptional('lat',       [],     @(X) isnumeric(X) && length(X) == 3);
    p.addParameter('loop',     '',     @(X) ischar(X) && any(strcmpi(X, {'component', 'none', ''})));
    p.addParameter('par',      false,  @isscalar);
    p.addParameter('output', []);
    p.addParameter('debug',  false);
    p.parse(ipsi, f, varargin{:});
    lat    = p.Results.lat;
    loop   = p.Results.loop;
    par    = p.Results.par;
    output = p.Results.output;
    debug  = p.Results.debug;
    
    if debug, fprintf('* pushImage\n'); end
    
    % --- Optimise parallelisation and splitting schemes
    [par, loop] = autoParLoop(par, loop, isa(f, 'file_array'), 1, size(f,4));
    
    % -- Read dim info
    dim = [size(ipsi) 1 1];
    nc  = size(f,4);
    
    % --- Default argument
    if isempty(lat)
        lat = dim(1:3);
    end
    
    % --- Allocate output
    if ~iscell(output)
        output = {output};
    end
    if numel(output) < 2
        output = [output {[]}];
    end
    
    ipsi = single(numeric(ipsi));
    
    % --- No loop
    if strcmpi(loop, 'none')
        if debug, fprintf('   - No loop\n'); end
        f = single(numeric(f));
        [pf1, c1] = spm_diffeo('pushc', f, ipsi, lat);
        if nargout > 2
            bb = boundingBox(c1);
        else
            bb.x = 1:size(c1, 1);
            bb.y = 1:size(c1, 2);
            bb.z = 1:size(c1, 3);
        end
        pf = prepareOnDisk(output{1}, [numel(bb.x) numel(bb.y) numel(bb.z) nc]);
        if nargout > 1
            c  = prepareOnDisk(output{2}, [numel(bb.x) numel(bb.y) numel(bb.z)]);
        end
        pf(:,:,:,:) = pf1(bb.x,bb.y,bb.z,:);
        clear pf1
        c(:,:,:)    = c1(bb.x,bb.y,bb.z);
        clear c1
        
        
    else % strcmpi(loop, 'component')
        if debug
            if par, fprintf('   - Parallelise on components\n');
            else,   fprintf('   - Serialise on components\n');   end
        end
        [pf1, c1] = spm_diffeo('pushc', single(f(:,:,:,1)), ipsi, lat);
        pf1(~isfinite(pf1)) = nan;
        c1(~isfinite(c1))   = nan;
        if nargout > 2
            bb = boundingBox(c1);
        else
            bb.x = 1:size(c1, 1);
            bb.y = 1:size(c1, 2);
            bb.z = 1:size(c1, 3);
        end
        pf = prepareOnDisk(output{1}, [numel(bb.x) numel(bb.y) numel(bb.z) nc]);
        if nargout > 1
            c  = prepareOnDisk(output{2}, [numel(bb.x) numel(bb.y) numel(bb.z)]);
        end
        pf(:,:,:,1) = pf1(bb.x,bb.y,bb.z);
        c(:,:,:)    = c1(bb.x,bb.y,bb.z);
        clear c1
        parfor (k=2:nc, par)
            pf1 = spm_diffeo('pushc', single(f(:,:,:,k)), ipsi, lat);
            pf1(~isfinite(pf1)) = nan;
            pf(:,:,:,k) = pf1(bb.x,bb.y,bb.z);
        end
        clear pf1
    end
        
    % --- Write on disk
    if ~isempty(output{1})
        pf = saveOnDisk(output{1}, pf, 'name', 'pf');
    end
    if nargout > 1 && ~isempty(output{2})
        c = saveOnDisk(output{2}, c, 'name', 'c');
    end
end
