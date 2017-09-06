classdef ImageSection < handle & AdvancedResize
    % (UI) Display a 2D section of an ND array.
    %% PUBLIC PROPERTIES
    properties (Dependent, AbortSet)
        ShowTicks
        FixedRatio
        Title
        ColorBar
        ColorBarPosition
        ColorMap
        ScaledColor
    end
    properties (AbortSet, SetObservable)
        Data
        VoxelSize
        WhichSection
        WhichPosition
        DefaultTitle
    end
    properties (SetAccess = private)
        OtherDims
    end
    %% PRIVATE PROPERTIES
    properties (Access = private)
        uiaxes
        uicolorbar
        uiimage
        uititle
        Section
        colorbarposition % To keep it stored when ther colorbar is hidden
        OldPosition
        OldSection
    end
    %% PUBLIC METHODS
    methods
        function obj = ImageSection(a, varargin)
            obj = obj@AdvancedResize(varargin);
            obj.Data = a;
            opt = getoptions(obj, varargin);
            obj.uipanel.Visible = 'off';
            
            obj.uiaxes = axes(...
                'Parent',               obj.uipanel, ...
                'Units',                'normalized', ...
                'OuterPosition',        [0 0 1 1]);
            
            % Properties needed to draw the image
            obj.WhichSection    = opt.WhichSection;
            obj.WhichPosition   = opt.WhichPosition;
            obj.VoxelSize       = opt.VoxelSize;
            obj.OtherDims       = opt.OtherDims;
            obj.DefaultTitle    = opt.DefaultTitle;
            
            extract_section(obj);
            if opt.ScaledColor
                obj.uiimage = image(obj.Section, 'Parent', obj.uiaxes, 'Visible', 'off', 'CDataMapping', 'Scaled');
            else
                obj.uiimage = image(obj.Section, 'Parent', obj.uiaxes, 'Visible', 'off');
            end
            if opt.FixedRatio
                daspect(obj.uiaxes, [obj.VoxelSize(obj.WhichSection) 1]);
            else
                daspect(obj.uiaxes, 'auto');
            end
            
            colormap(obj.uiaxes, opt.ColorMap);
            obj.uicolorbar = 0;
            if opt.ColorBar
                obj.uicolorbar = colorbar(obj.uiaxes, ...
                    'Location', opt.ColorBarPosition);
            end
            obj.colorbarposition = opt.ColorBarPosition;
            obj.uititle = 0;
            if ~isempty(opt.Title)
                obj.uititle = title(obj.uiaxes, opt.Title);
            end
            
            addlistener(obj, 'Data', 'PostSet', @(s,e)updatedata(obj));
            addlistener(obj, 'VoxelSize', 'PostSet', @(s,e)updatevs(obj));
            addlistener(obj, 'WhichSection', 'PreSet', @(s,e)updatesectionpre(obj));
            addlistener(obj, 'WhichSection', 'PostSet', @(s,e)updatesectionpost(obj));
            addlistener(obj, 'WhichPosition', 'PostSet', @(s,e)updateposition(obj));

            obj.uipanel.Visible = 'on';
            obj.uiimage.Visible = 'on';
        end
        %% SETTERS/GETTERS
        function value = get.ShowTicks(obj)
            value = ~isempty(obj.uiaxes.XTickLabel);
        end
        function set.ShowTicks(obj, value)
            if value
                obj.uiaxes.XTickLabelMode = 'auto';
                obj.uiaxes.YTickLabelMode = 'auto';
            else
                obj.uiaxes.XTick = [];
                obj.uiaxes.XTickLabel = [];
                obj.uiaxes.YTick = [];
                obj.uiaxes.YTickLabel = [];
            end
            drawnow
        end
        function value = get.FixedRatio(obj)
            value = strcmp(daspect(obj.uiaxes, 'mode'), 'manual');
        end
        function set.FixedRatio(obj, value)
            if value
                daspect(obj.uiaxes, [obj.VoxelSize(obj.WhichSection) 1]);
            else
                daspect(obj.uiaxes, 'auto');
            end
            drawnow
        end
        function value = get.Title(obj)
            if obj.uititle ~= 0
                value = obj.uititle.String;
            else
                value = '';
            end
        end
        function set.Title(obj, value)
            autotitle(obj, value);
        end
        function value = get.ColorBar(obj)
            value = obj.uicolorbar ~= 0;
        end
        function set.ColorBar(obj, value)
            colorbar(obj.uiaxes, 'off');
            if value
                obj.uicolorbar = colorbar(obj.uiaxes, ...
                    'Location', obj.colorbarposition);
            else
                obj.uicolorbar = 0;
            end
            drawnow
        end
        function value = get.ColorBarPosition(obj)
            if obj.uicolorbar ~= 0
                value = obj.uicolorbar.Location;
            else
                value = obj.colorbarposition;
            end
        end
        function set.ColorBarPosition(obj, value)
            if obj.uicolorbar ~= 0
                obj.uicolorbar.Location = value;
            end
            obj.colorbarposition = value;
            drawnow
        end
        function value = get.ColorMap(obj)
            value = obj.uicolormap;
        end
        function set.ColorMap(obj, value)
            obj.uicolormap = colormap(obj.uiaxes, value);
            drawnow
        end
        function value = get.ScaledColor(obj)
            value = obj.uiimage.CDataMapping == 'scaled';
        end
        function set.ScaledColor(obj, value)
            if value
                obj.uiimage.CDataMapping = 'scaled';
            else
                obj.uiimage.CDataMapping = 'direct';
            end
            drawnow
        end
    end
    %% PRIVATE METHODS
    methods(Access = private)
        %% UTILS
        function opt = getoptions(obj, in)
            % Default values
            opt                     = struct;
            opt.ShowTicks           = true;
            opt.FixedRatio          = true;
            opt.Title               = '';
            opt.DefaultTitle        = true;
            opt.ColorBar            = true;
            opt.ColorBarPosition    = 'eastoutside';
            opt.ColorMap            = 'default';
            opt.ScaledColor         = true;
            opt.WhichSection        = [1 2];
            opt.VoxelSize           = ones(size(obj.Data));
            
            % Parse
            opt = parse_varargin(in, opt);
            
            % Post processing
            if opt.DefaultTitle && isempty(opt.Title)
                if issame(opt.WhichSection, [1 2])
                    opt.Title = 'Coronal';
                elseif issame(opt.WhichSection, [1 3])
                    opt.Title = 'Axial';
                elseif issame(opt.WhichSection, [2 3])
                    opt.Title = 'Sagittal';
                end
            else
                opt.Title = char(opt.Title);
            end
            
            dim = size(obj.Data);
            opt.OtherDims = [];
            if ~isfield(opt, 'WhichPosition')
                opt.WhichPosition = [];
                do_position = true;
            else
                do_position = false;
            end
            for i=1:length(dim)
                if ~numel(find(opt.WhichSection == i))
                    opt.OtherDims = [opt.OtherDims i];
                    if do_position
                        opt.WhichPosition = [opt.WhichPosition ceil(dim(i)/2)];
                    end
                end
            end
        end
        %%
        function extract_section(obj)
            obj.Section = obj.Data;
            dim = size(obj.Data);
            for k=1:length(obj.OtherDims)
                obj.Section = arrayutils('select_slice', obj.Section, ...
                    obj.WhichPosition(k), obj.OtherDims(k));
            end
            % I am doing a cast towards single as a (temporar) way to make
            % it work with handle_array objects.
            obj.Section = reshape(numeric(obj.Section), dim(obj.WhichSection));
        end
        function autotitle(obj, value)
            if nargin < 2
                value = obj.Title;
            elseif isempty(value) && obj.uititle ~= 0
                delete(obj.uititle);
                obj.uititle = 0;
            end
            if obj.DefaultTitle && isempty(value)
                if issame(obj.WhichSection, [1 2])
                    newtitle = 'Coronal';
                elseif issame(obj.WhichSection, [1 3])
                    newtitle = 'Axial';
                elseif issame(obj.WhichSection, [2 3])
                    newtitle = 'Sagittal';
                else
                    newtitle = ['[', ...
                                char(obj.WhichSection(1)), ...
                                ' ', ...
                                char(obj.WhichSection(2)), ... 
                                ']'];
                end
            else
                newtitle = char(value);
            end
            obj.uititle = title(obj.uiaxes, newtitle);
            drawnow
        end
        %% CALLBACKS
        function updatedata(obj)
            extract_section(obj);
            obj.uiimage.CData = obj.Section;
            drawnow
        end
        function updatevs(obj)
            if obj.FixedRatio
                daspect(obj.uiaxes, [obj.VoxelSize(obj.WhichSection) 1]);
            else
                daspect(obj.uiaxes, 'auto');
            end
            drawnow
        end
        function updatesectionpre(obj)
            obj.OldPosition = obj.WhichPosition;
            obj.OldSection  = obj.WhichSection;
        end
        function updatesectionpost(obj)
            dim = size(obj.Data);
            otherDims = zeros([1 length(dim)-2]);
            whichPosition = zeros([1 length(dim)-2]);
            k = 1;
            l = 1;
            for i=1:length(dim)
                if ~numel(find(obj.WhichSection == i))
                    otherDims(k) = i;
                    if ~numel(find(obj.OldSection == i))
                        whichPosition(k) = obj.OldPosition(l);
                        l = l+1;
                    else
                        whichPosition(k) = ceil(dim(i)/2);
                    end
                    k = k+1;
                elseif ~numel(find(obj.OldSection == i))
                    l = l+1;
                end
            end
            obj.OtherDims = otherDims;
            obj.WhichPosition = whichPosition;
            updatedata(obj);
            updatevs(obj);
            autotitle(obj);
            drawnow
        end
        function updateposition(obj)
            updatedata(obj);
            updatevs(obj);
            autotitle(obj);
            drawnow
        end
    end
end