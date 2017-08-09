classdef SectionSlider < AdvancedResize & handle
    
    %% PUBLIC PROPERTIES
    properties (Dependent, AbortSet)
        Coordinates
        Dims
        Length
    end
    properties (SetAccess = private, SetObservable)
        CoordinatesChanged
    end
    properties (Dependent, AbortSet)
        SliderHeight
    end
    
    %% PRIVATE PROPERTIES
    properties (Access = public)
        sliders
    end
    
    %% PUBLIC METHODS
    methods
        %% CONSTURCTOR
        function obj = SectionSlider(varargin)
            obj = obj@AdvancedResize(varargin);
            obj.uipanel.BackgroundColor = 'red';
            opt = getoptions(obj, varargin);
            obj.uipanel.Visible = 'off';
            obj.HeightUnit = opt.HeightUnit;
            obj.Height = opt.SliderHeight * length(opt.Dims);
            
            % Callbacks
            obj.CoordinatesChanged = false;
            
            obj.sliders = cell(length(opt.Dims), 1);
            for i=1:length(opt.Dims)
                obj.sliders{i} = DiscreteSlider( ...
                    'Parent',               obj.uipanel, ...
                    'HeightUnit',           obj.HeightUnit, ...
                    'Height',               opt.SliderHeight, ...
                    'Length',               opt.Length(i), ...
                    'Name',                 opt.Dims(i), ...
                    'VerticalOrigin',       (i-1)*opt.SliderHeight, ...
                    'VerticalOriginSide',   'top', ...
                    'VerticalOriginUnit',   obj.HeightUnit, ...
                    'Visible',              'off');
                if isfield(opt, 'Coordinates')
                    obj.sliders{i}.Value = opt.Coordinates(i);
                end
                addlistener(obj.sliders{i}, 'ValueChanged', 'PostSet', @(s,e)changeCoordinates(obj));
            end
            
            
            % Draw
            obj.uipanel.Visible = 'on';
            for i=1:length(opt.Dims)
                obj.sliders{i}.Visible = 'on';
            end
        end
        %% GETTERS/SETTERS
        function value = get.Coordinates(obj)
            value = zeros([1 length(obj.sliders)]);
            for i=1:length(obj.sliders)
                value(i) = obj.sliders{i}.Value;
            end 
        end
        function set.Coordinates(obj, value)
            for i=1:length(obj.sliders)
                obj.sliders{i}.Value = value(i);
            end 
            drawnow
        end
        function value = get.Dims(obj)
            value = zeros([1 length(obj.sliders)]);
            for i=1:length(obj.sliders)
                v = obj.sliders{i}.Name;
                switch v
                    case 'X'
                        v = 1;
                    case 'Y'
                        v = 2;
                    case 'Z'
                        v = 3;
                    otherwise
                        v = str2double(v);
                end
                value(i) = v;
            end 
        end
        function set.Dims(obj, value)
            for i=1:length(obj.sliders)
                obj.sliders{i}.Name = value(i);
            end 
            drawnow
        end
        function value = get.Length(obj)
            value = zeros([1 length(obj.sliders)]);
            for i=1:length(obj.sliders)
                value(i) = obj.sliders{i}.Length;
            end 
        end
        function set.Length(obj, value)
            for i=1:length(obj.sliders)
                obj.sliders{i}.Length = value(i);
            end 
            drawnow
        end
        function value = get.SliderHeight(obj)
            value = obj.sliders{1}.Height;
        end
        function set.SliderHeight(obj, value)
            for i=1:length(obj.sliders)
                obj.sliders{i}.Height = value;
                obj.sliders{i}.VerticalOrigin = (i-1)*SliderHeigh;
            end 
            obj.Height = length(obj.sliders) * value;
            drawnow
        end
    end
    
    %% PRIVATE METHODS
    methods (Access = private)
        function opt = getoptions(~, in)
            % Default values
            opt.Dims = 1;
            opt.Length = 100;
            opt.SliderHeight = 20;
            opt.HeightUnit = 'pixels';
            
            % Parse
            opt = parse_varargin(in, opt);
        end
        function changeCoordinates(obj)
            obj.CoordinatesChanged = ~obj.CoordinatesChanged;
        end
    end
end