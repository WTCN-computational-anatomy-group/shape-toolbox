function animateDiffeo(iphi, varargin)
% _________________________________________________________________________
%
% FORMAT animateDiffeo(iphi, (a), (iscat), ...)
%
% -------------------------------------------------------------------------
%
% DATA
% ----
% iphi   - Diffeomorphic flow (obtained from shootGeodesicInv)
%
% OPTIONAL
% --------
% a      - Image to deform [none]
% iscat  - Image to deform is categorical and should be softmaxed [false]
%
% KEYWORD ARGUMENTS
% -----------------
% def        - Deformation video (output path) [do not generate]
% grid       - Grid deformation video (output path) [do not generate]
% image      - Image deformation video (output path) [do not generate]
% slice      - Slice direction (1 = sagittal, 2 = coronal, [3] = axial)
% x          - Slice location [sagittal: 1/3, axial/coronal: middle]
% vs         - Lattice voxel size [1 1 1]
% grid_space - Nb of voxels between lines in the grid [10 10 10]
% loop       - Repeat frames to form a loop? [false]
% duration   - Total duration of the video in sec [5]
%
% Note: Possible output extensions are .mp4, .avi, .gif 
% _________________________________________________________________________

    % Parse arguments
    % ---------------
    p = inputParser;
    p.FunctionName = 'animateDiffeo';
    p.addRequired('iphi',   @(X) isnumeric(X) && size(X,5) > 1);
    p.addOptional('a',      []);
    p.addOptional('iscat',  false, @(X) isscalar(X) && (isnumeric(X)||islogical(X)));
    p.addParameter('def',   '',    @ischar);
    p.addParameter('grid',  '',    @ischar);
    p.addParameter('image', '',    @ischar);
    p.addParameter('slice', 3,     @(X) isscalar(X) && isnumeric(X));
    p.addParameter('x',     nan,   @(X) isscalar(X) && isnumeric(X));
    p.addParameter('lat',   [],    @isnumeric);
    p.addParameter('vs',  [1 1 1], @isnumeric);
    p.addParameter('grid_space', [10 10 10], @isnumeric);
    p.addParameter('loop',  false, @(X) isscalar(X) && (isnumeric(X)||islogical(X)));
    p.addParameter('duration', 5,     @(X) isscalar(X) && isnumeric(X));
    p.parse(iphi, varargin{:});
    a     = p.Results.a;
    iscat = p.Results.iscat;
    def   = p.Results.def;
    grid  = p.Results.grid;
    image = p.Results.image;
    slice = p.Results.slice;
    x     = p.Results.x;
    vs    = p.Results.vs;
    grid_space = p.Results.grid_space;
    loop       = p.Results.loop;
    duration   = p.Results.duration;
    
    dim = [size(iphi) 1 1 1];
    lat = dim(1:3);
    if ~isempty(grid)
        grid_image = zeros(lat, 'single');
        grid_image(grid_space(1):grid_space(1):end,:,:) = 1;
        grid_image(:,grid_space(2):grid_space(2):end,:) = 1;
        grid_image(:,:,grid_space(3):grid_space(3):end) = 1;
    end
    
    if ~isfinite(x)
        switch slice
            case {2,3}
                x = ceil(dim(slice)/2);
            case 1
                x = ceil(dim(slice)/3);
        end
    end

    N   = size(iphi, 5);
    fig = figure;
    if ~isempty(def)
        frames_def(1:N) = struct('cdata',[],'colormap',[]);
        id = spm_warps('identity', lat);
        nrm = max(max(max(max(max(sum((iphi-id).^2, 4))))));
    end
    if ~isempty(grid)
        frames_grid(1:N)  = struct('cdata',[],'colormap',[]);
    end
    if ~isempty(image)
        frames_image(1:N) = struct('cdata',[],'colormap',[]);
    end
    fprintf('Build and record frames\n');
    for i=1:N
        fprintf('i = %3d', i);
        
        if ~isempty(def)
            fprintf(' | Def');
            figure(fig);
            h = imshow_deformation(iphi(:,:,:,:,i)-id, vs, x, slice, nrm);
            frames_def(i) = getframe(h.Parent);
        end
        
        if ~isempty(grid)
            fprintf(' | Grid');
%             g1 = warp(iphi(:,:,:,:,i), grid_image);
%             switch slice
%                 case 1
%                     g1 = reshape(g1(x,:,:), [dim(2) dim(3)]);
%                 case 2
%                     g1 = reshape(g1(:,x,:), [dim(1) dim(3)]);
%                 case 3
%                     g1 = reshape(g1(:,:,x), [dim(1) dim(2)]);
%                 otherwise
%                     error('not handled')
%             end
%             g1 = permute(g1, [2 1]);
%             figure(fig);
%             h = imagesc(g1);
%             colormap(gray)
            figure(fig);
            switch slice
                case 1
                    g11 = reshape(iphi(x,:,:,2,i), [dim(2) dim(3)]);
                    g12 = reshape(iphi(x,:,:,3,i), [dim(2) dim(3)]);
                    xmin = min(min(min(iphi(x,:,:,2,:))));
                    xmax = max(max(max(iphi(x,:,:,2,:))));
                    ymin = min(min(min(iphi(x,:,:,3,:))));
                    ymax = max(max(max(iphi(x,:,:,3,:))));
                case 2
                    g11 = reshape(iphi(:,x,:,1,i), [dim(1) dim(3)]);
                    g12 = reshape(iphi(:,x,:,3,i), [dim(1) dim(3)]);
                    xmin = min(min(min(iphi(:,x,:,1,:))));
                    xmax = max(max(max(iphi(:,x,:,1,:))));
                    ymin = min(min(min(iphi(:,x,:,3,:))));
                    ymax = max(max(max(iphi(:,x,:,3,:))));
                case 3
                    g11 = reshape(iphi(:,:,x,1,i), [dim(1) dim(2)]);
                    g12 = reshape(iphi(:,:,x,2,i), [dim(1) dim(2)]);
                    xmin = min(min(min(iphi(:,:,x,1,:))));
                    xmax = max(max(max(iphi(:,:,x,1,:))));
                    ymin = min(min(min(iphi(:,:,x,2,:))));
                    ymax = max(max(max(iphi(:,:,x,2,:))));
                otherwise
                    error('not handled')
            end
            g11 = g11(1:5:end,1:5:end);
            g12 = g12(1:5:end,1:5:end);
%             g11 = permute(g11, [2 1]);
%             g12 = permute(g12, [2 1]);
            h = plot(g11,g12,'k',g11',g12','k');
            xlim([xmin xmax])
            ylim([ymin ymax])
            axis off
            switch slice
                case 1
                    daspect([1/vs(2) 1/vs(3) 1]);
                case 2
                    daspect([1/vs(1) 1/vs(3) 1]);
                case 3
                    daspect([1/vs(1) 1/vs(2) 1]);
            end
            frames_grid(i) = getframe(h(1).Parent);
            clear g1
        end
            
        if ~isempty(image)
            fprintf(' | Image');
            a1 = warp(iphi(:,:,:,:,i), a);
            if iscat
                a1 = reconstructProbaTemplate(a1);
            end
            figure(fig);
            h = imshow_categorical(a1, vs, x, slice);
            frames_image(i) = getframe(h.Parent);
            clear a1
        end
        
        fprintf('\n');
        
    end
    close(fig);
    
    fprintf('Save frames\n');
    
    for j=1:3
        switch j
            case 1
                fname  = def;
                if isempty(fname)
                    continue
                end
                frames = frames_def;
            case 2
                fname  = grid;
                if isempty(fname)
                    continue
                end
                frames = frames_grid;
            case 3
                fname  = image;
                if isempty(fname)
                    continue
                end
                frames = frames_image;
        end
        fprintf('%s\n', fname);
    
        if loop
            iter = [1:N N-1:-1:2];
        else
            iter = 1:N;
        end
        framerate = numel(iter)/duration;
        if endsWith(fname, 'mp4'), vidopt = {'MPEG-4'};
        else,                      vidopt = {};
        end
        if endsWith(fname, 'avi') || endsWith(fname, 'mp4')
            vid = VideoWriter(fname, vidopt{:});
            vid.FrameRate = framerate;
            open(vid);
        end
        first_frame = true;
        for i=iter
            if i > 1
                dim1 = size(frames(i).cdata);
                dim0 = size(frames(1).cdata);
                if dim1(1) > dim0(1)
                    frames(i).cdata = frames(i).cdata(1:dim0(1),:,:);
                end
                if dim1(2) > dim0(2)
                    frames(i).cdata = frames(i).cdata(:,1:dim0(2),:);
                end
                if dim1(1) < dim0(1)
                    frames(i).cdata = padarray(frames(i).cdata, [dim0(1)-dim1(1) 0 0], 0, 'post');
                end
                if dim1(2) < dim0(2)
                    frames(i).cdata = padarray(frames(i).cdata, [0 dim0(2)-dim1(2) 0], 0, 'post');
                end
            end
            if endsWith(fname, 'gif')
                [imind, cm] = rgb2ind(frame2im(frames(i)), 256);
                if first_frame
                    imwrite(imind, cm, fname, 'gif', 'Loopcount', inf, 'DelayTime', 1/framerate); 
                    first_frame = false;
                else
                    imwrite(imind, cm, fname, 'gif', 'WriteMode', 'append', 'DelayTime', 1/framerate); 
                end
            elseif endsWith(fname, 'avi') || endsWith(fname, 'mp4') 
                writeVideo(vid, frames(i));
            end
        end
        if endsWith(fname, 'avi') || endsWith(fname, 'mp4') 
            close(vid);
        end
    end

end