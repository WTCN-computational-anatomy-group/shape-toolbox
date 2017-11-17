classdef DiscreteSlider < handle & AdvancedResize
% (UI) Discretized slider with additional properties
    %% PUBLIC PROPERTIES
    properties (Dependent, AbortSet)
        ShowName        % Show the slider name
        ShowValue       % Show the slider position in a box
        EditableValue   % Is the slider position editable through the box
        Name
    end
    properties (Dependent, SetAccess = public)
        Length          % Alternative way of setting max - default: 100
    end
    properties (Dependent, AbortSet, SetObservable)
        Min             % default: 1
        Max             % default: length
        Value           % default: closest(min + (max-min)/2)
    end
    properties (SetAccess = private, SetObservable)
        MinChanged
        MaxChanged
        ValueChanged
    end
    properties (AbortSet, SetObservable)
        Step            % default: 1
    end
    properties
        DefaultName     % default: true
    end
    %% PRIVATE PROPERTIES
    properties (Access = public)
        uilabel     % Name of the slider
        uislider    % Slider
        uivalue     % Editable value
    end
    %% PUBLIC METHODS
    methods
        %% CONSTRUCTOR
        function obj = DiscreteSlider(varargin)
            obj = obj@AdvancedResize(varargin);
            opt = getoptions(obj, varargin);
            obj.uipanel.Visible = 'off';
            
            obj.uilabel = uicontrol( ...
                'Parent',           obj.uipanel, ...
                'Style',            'text', ...
                'Unit',             'pixels', ...
                'String',           opt.Name, ...
                'Visible',          'off');
            obj.uislider = uicontrol( ...
                'Parent',           obj.uipanel, ...
                'Style',            'slider', ...
                'Unit',             'pixels', ...
                'Value',            opt.Value, ...
                'Min',              opt.Min, ...
                'Max',              opt.Max, ...
                'Visible',          'off');
            obj.uivalue = uicontrol( ...
                'Parent',           obj.uipanel, ...
                'Style',            'edit', ...
                'Unit',             'pixels', ...
                'Visible',          'off');
            obj.DefaultName     = opt.DefaultName;   % Non-dependent property
            obj.Step            = opt.Step;          % Non-dependent property
            obj.EditableValue   = opt.EditableValue; % Triggers ui change
            obj.Length          = opt.Length;        % Triggers ui change
            % Listeners
            addlistener(obj, 'Step', 'PostSet', @(s,e)updatestep(obj));
            addlistener(obj, 'Min', 'PostSet', @(s,e)updatestep(obj));
            addlistener(obj, 'Max', 'PostSet', @(s,e)updatestep(obj));
            % Callbacks
            obj.MinChanged = false;
            obj.MaxChanged = false;
            obj.ValueChanged = false;
            addCallback(obj, @(s,e)resizecontent(obj));
            obj.uislider.Callback = @(s,e)sliderchanged(obj);
            obj.uivalue.Callback = @(s,e)valuechanged(obj);
            % Draw
            updatestep(obj);
            resizecontent(obj);
            if opt.ShowName
                obj.uilabel.Visible = 'on';
            else
                obj.uilabel.Visible = 'off';
            end
            obj.uislider.Visible = 'on';
            if opt.ShowValue
                obj.uivalue.Visible = 'on';
            else
                obj.uivalue.Visible = 'off';
            end
            sliderchanged(obj);
            obj.uipanel.Visible = 'on';
        end
        %% SETTERS/GETTERS
        function set.Min(obj, value)
            obj.uislider.Min = value;
            obj.MinChanged = ~obj.MinChanged;
            drawnow
        end
        function value = get.Min(obj)
            value = obj.uislider.Min;
        end
        function set.Max(obj, value)
            obj.uislider.Max = value;
            obj.MaxChanged = ~obj.MaxChanged;
            drawnow
        end
        function value = get.Max(obj)
            value = obj.uislider.Max;
        end
        function set.Length(obj, value)
            obj.uislider.Max = obj.uislider.Min + value - obj.Step;
            drawnow
        end
        function set.Value(obj, value)
            obj.uislider.Value = value;
            obj.ValueChanged = ~obj.ValueChanged;
            drawnow
        end
        function value = get.Value(obj)
            value = obj.uislider.Value;
        end
        function set.Name(obj, value)
            if obj.DefaultName
                switch value
                    case 1
                        value = 'X';
                    case 2
                        value = 'Y';
                    case 3
                        value = 'Z';
                    otherwise
                        value = char(value);
                end
            else
                value = char(value);
            end
            obj.uilabel.String = value;
        end
        function value = get.Name(obj)
            value = obj.uilabel.String;
        end
        function set.ShowName(obj, value)
            if value
                obj.uilabel.Visible = 'on';
            else
                obj.uilabel.Visible = 'off';
            end
            drawnow
        end
        function value = get.ShowName(obj)
            value = strcmp(obj.uilabel.Visible, 'on');
        end
        function set.ShowValue(obj, value)
            if value
                obj.uivalue.Visible = 'on';
            else
                obj.uivalue.Visible = 'off';
            end
            drawnow
        end
        function value = get.ShowValue(obj)
            value = strcmp(obj.uivalue.Visible, 'on');
        end
        function set.EditableValue(obj, value)
            if value
                obj.uivalue.Style = 'edit';
            else
                obj.uivalue.Style = 'text';
            end
            drawnow
        end
        function value = get.EditableValue(obj)
            value = obj.uivalue.Style == 'edit';
        end
    end
    %% PRIVATE METHODS
    methods (Access = private)
        %% CALLBACKS
        function sliderchanged(obj)
            value = obj.uislider.Value;
            value = obj.Min + obj.Step * round((value - obj.Min) / obj.Step);
            if value < obj.Min
                value = obj.Min;
            elseif value > obj.Max
                value = obj.Max;
            end
            obj.uislider.Value = value;
            obj.uivalue.String = num2str(value);
            obj.ValueChanged = ~obj.ValueChanged;
        end
        function valuechanged(obj)
            value = str2double(obj.uivalue.String);
            value = obj.Min + obj.Step * round((value - obj.Min) / obj.Step);
            if value < obj.Min
                value = obj.Min;
            elseif value > obj.Max
                value = obj.Max;
            end
            obj.uivalue.String = num2str(value);
            obj.uislider.Value = value;
            obj.ValueChanged = ~obj.ValueChanged;
        end
        function resizecontent(obj)
            % panel
            panel_units = obj.uipanel.Units;
            obj.uipanel.Units = 'pixels';
            panel_pos = obj.uipanel.Position;
            % fixed sizes
            if obj.ShowName
                width_name = 25;
            else
                width_name = 0;
            end
            if obj.ShowValue
                width_value = 25;
            else
                width_value = 0;
            end
            % label
            label_units = obj.uilabel.Units;
            obj.uilabel.Units = 'pixels';
            obj.uilabel.Position = [0 ...
                                    0 ...
                                    width_name ...
                                    panel_pos(4)];
            % slider
            slider_units = obj.uislider.Units;
            obj.uislider.Units = 'pixels';
            obj.uislider.Position = [obj.uilabel.Position(1) + obj.uilabel.Position(3) ...
                                     0 ...
                                     panel_pos(3) - width_name - width_value ...
                                     panel_pos(4)];
            % slider
            value_units = obj.uivalue.Units;
            obj.uivalue.Units = 'pixels';
            obj.uivalue.Position = [obj.uislider.Position(1) + obj.uislider.Position(3) ...
                                     0 ...
                                     width_value ...
                                     panel_pos(4)];
            %% end
            obj.uilabel.Units = label_units;
            obj.uislider.Units = slider_units;
            obj.uivalue.Units = value_units;
            obj.uipanel.Units = panel_units;
        end
        function updatestep(obj)
            if obj.Max == obj.Min
                step = 1;
            else
                step = obj.Step/(obj.Max - obj.Min);
            end
            obj.uislider.SliderStep = [step max(0.1, step)];
            drawnow
        end
    end
    methods (Access = private, Static)
        %% UTILS
        function opt = getoptions(~, in)
            % Default values
            opt                     = struct;
            opt.Min                 = 1;
            opt.Max                 = 100;
            opt.Step                = 1;
            opt.Name                = '';
            opt.DefaultName         = true;
            opt.ShowName            = 'on';
            opt.ShowValue           = 'on';
            opt.EditableValue       = true;
            
            % Parse
            opt = parse_varargin(in, opt);
            
            % Post processing
            if opt.DefaultName
                switch opt.Name
                    case 1
                        opt.Name = 'X';
                    case 2
                        opt.Name = 'Y';
                    case 3
                        opt.Name = 'Z';
                    otherwise
                        opt.Name = char(opt.Name);
                end
            else
                opt.Name = char(opt.Name);
            end
            if isfield(opt, 'Length')
                opt.Max = opt.Min + opt.Length - 1;
            else
                opt.Length = opt.Max - opt.Min + 1;
            end
            if ~isfield(opt, 'Value')
                opt.Value = opt.Min + opt.Step * round((opt.Max - opt.Min) / 2 / opt.Step);
            end
        end
    end
end