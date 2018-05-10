function varargout = shape_ui(id, varargin)

    switch lower(id)
        case 'title'
            [varargout{1:nargout}] = title(varargin{:});
        case 'posttitle'
            [varargout{1:nargout}] = posttitle(varargin{:});
        case 'em'
            [varargout{1:nargout}] = em(varargin{:});
        case 'lb'
            [varargout{1:nargout}] = lb(varargin{:});
        case 'postlb'
            [varargout{1:nargout}] = postlb(varargin{:});
    end

end

function start = title(text, newline, time)
% FORMAT start = shape_ui('Title', text, newline, time)
    if nargin < 3
        time = false;
        if nargin < 2
            newline = true;
            if nargin < 1
                text = '';
            end
        end
    end
    if time
        start = tic;
    else
        start = 0;
    end
    fprintf('%-27s | ', text);
    if newline
        fprintf('\n');
    end
end

function posttitle(time)
% FORMAT shape_ui('PostTitle', time)
    fprintf('%s\n', sec2ydhms(time, true));
end

function em(iteration, text)
% FORMAT shape_ui('EM', iteration, text)

    if nargin < 2
        text = 'VEM';
        if nargin < 1
            iteration = 'Init';
        end
    end
    
    if ischar(iteration)
        fprintf('%-20s%7s | %s |\n', text, iteration, repmat('=',1,48));
    else
        fprintf('%-20s%7d | %s |\n', text, iteration, repmat('=',1,48));
    end
end

function val0 = lb(val0, text)
% FORMAT val0 = shape_ui('LB', val, text)

    if nargin < 2
        text = 'Update Lower Bound';
    end
    fprintf('%-27s | %10.3g -> ', text, val0);
end


function postlb(val, val0)
% FORMAT val0 = shape_ui('PostLB', val, val0)

    if val > val0
        sign = '(+)';
    elseif val < val0
        sign = '(-)';
    else
        sign = '(=)';
    end
    fprintf('%10.3g %s\n', val, sign);
end