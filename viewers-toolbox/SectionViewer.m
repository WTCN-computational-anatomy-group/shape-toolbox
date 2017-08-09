classdef SectionViewer < AdvancedResize & handle
    %% PUBLIC PROPERTIES
    properties (Dependent = true, AbortSet = true)
        WhichPosition
        WhichSection
        Title
    end
    %% PRIVATE PROPERTIES
    properties (Access = public) % private
        imageSection
        sectionSlider
    end
    %% PUBLIC METHODS
    methods
        %% CONSTRUCTOR
        function obj = SectionViewer(a, varargin)
            obj = obj@AdvancedResize(varargin);
            opt = getoptions(obj, a, varargin);
            obj.uipanel.Visible = 'off';
            
            obj.sectionSlider = SectionSlider( ...
                'Parent',               obj.uipanel, ...
                'Dims',                 opt.Dims, ...
                'Length',               opt.Length);
            % init position
            otherdims = [];
            for i=1:length(size(a))
                if ~numel(find(opt.WhichSection == i))
                    otherdims = [otherdims i];
                end
            end
            pos = SectionViewer.slider2ImageCoord(...
                    obj.sectionSlider.Coordinates, size(a), otherdims); 
            obj.imageSection = ImageSection(a, ...
                'Parent',               obj.uipanel, ...
                'WhichSection',         opt.WhichSection, ...
                'WhichPosition',        pos, ...
                'VerticalOriginSide',   'top', ...
                'HeightOffset',         obj.sectionSlider.Height);
            
            addlistener(obj.sectionSlider, 'CoordinatesChanged', 'PostSet', @(s,e) slider2image(obj));
            if isfield(opt, 'WhichPosition')
                obj.imageSection.WhichPosition = opt.WhichPosition;
            end
            if isfield(opt, 'Title')
                obj.imageSection.Title = opt.Title;
            end
                
            obj.uipanel.Visible = 'on';
        end
        %% GETTERS/SETTERS
        function value = get.WhichPosition(obj)
            value = obj.imageSection.WhichPosition;
        end
        function set.WhichPosition(obj, value)
            obj.imageSection.WhichPosition = value;
        end
        function value = get.WhichSection(obj)
            value = obj.imageSection.WhichSection;
        end
        function set.WhichSection(obj, value)
            obj.imageSection.WhichSection = value;
        end
        function value = get.Title(obj)
            value = obj.imageSection.Title;
        end
        function set.Title(obj, value)
            obj.imageSection.Title = value;
        end
    end
    %% PRIVATE METHODS
    methods (Access = private)
        %% UTILS
        function opt = getoptions(~, a, in)
            % Default values
            dim = size(a);
            opt                     = struct;
            opt.WhichSection        = [1 2];
            opt.VoxelSize           = ones([1 length(dim)]);
            
            % Parse
            opt = parse_varargin(in, opt);
            
            opt.Length = [];
            opt.Dims = [];
            opt.Names = [];
            for i=1:length(dim)
                if ~numel(find(opt.WhichSection == i)) && dim(i) > 1
                    opt.Dims = [opt.Dims i];
                    opt.Length = [opt.Length dim(i)];
                end
            end
        end
        %% CALLBACKS
        function slider2image(obj)
            pos = SectionViewer.slider2ImageCoord(...
                    obj.sectionSlider.Coordinates, ...
                    size(obj.imageSection.Data), ...
                    obj.imageSection.OtherDims);
            obj.imageSection.WhichPosition = pos;
        end
            
    end
    %% STATIC HELPER
    methods (Access = private, Static)
        function pos = slider2ImageCoord(coord, alldims, otherdims)
            pos = ones(size(otherdims));
            dims = alldims(otherdims);
            l = 1;
            for i=1:length(otherdims)
                if dims(i) > 1
                    pos(i) = coord(l);
                    l = l + 1;
                end
            end
        end
    end
end