function varargout = dt(id, varargin)
% FORMAT Y = dt('fft', X, d)
% FORMAT Y = dt('fftn', X)
% FORMAT Y = dt('ifft', X, d)
% FORMAT Y = dt('ifftn', X)
% FORMAT Y = dt('dct', X, d)
% FORMAT Y = dt('dctn', X)
% FORMAT Y = dt('idct', X, d)
% FORMAT Y = dt('idctn', X)
% FORMAT Y = dt('dst', X, d)
% FORMAT Y = dt('dstn', X)
% FORMAT Y = dt('idst', X, d)
% FORMAT Y = dt('idstn', X)
%
% Collection of tools to compute discrete spatial/frequency transforms.
%
% FORMAT help dt>function_
% Returns the help file of the selected function.
    switch id
        case 'fft'
            [varargout{1:nargout}] = fft_(varargin{:});
        case 'fftn'
            [varargout{1:nargout}] = fftn_(varargin{:});
        case 'ifft'
            [varargout{1:nargout}] = ifft_(varargin{:});
        case 'ifftn'
            [varargout{1:nargout}] = ifftn_(varargin{:});
        case 'dct'
            [varargout{1:nargout}] = dct_(varargin{:});
        case 'dctn'
            [varargout{1:nargout}] = dctn_(varargin{:});
        case 'idct'
            [varargout{1:nargout}] = idct_(varargin{:});
        case 'idctn'
            [varargout{1:nargout}] = idctn_(varargin{:});
        case 'dst'
            [varargout{1:nargout}] = dst_(varargin{:});
        case 'dstn'
            [varargout{1:nargout}] = dstn_(varargin{:});
        case 'idst'
            [varargout{1:nargout}] = idst_(varargin{:});
        case 'idstn'
            [varargout{1:nargout}] = idstn_(varargin{:});
        otherwise
            error('No function named %s\n', string(id));
    end

% === Subfunctions ========================================================

function X = fft_(X, d)
    if nargin < 2
        X = fft(X);
    else
        X = fft(X, [], d);
    end
    
function X = fftn_(X)
    X = fftn(X);

function X = ifft_(X, d)
    if nargin < 2
        X = ifft(X);
    else
        X = ifft(X, [], d);
    end
    
function X = ifftn_(X)
    X = ifftn(X);
    
function X = dct_(X, d)
    if nargin < 2
        d = 1;
    end
    X = make_d(@dct, X, d);
    
function X = dctn_(X)
    X = make_n(@dct, X);

function X = idct_(X, d)
    if nargin < 2
        d = 1;
    end
    X = make_d(@idct, X, d);
    
function X = idctn_(X)
    X = make_n(@idct, X);
    
function X = dst_(X, d)
    if nargin < 2
        d = 1;
    end
    X = make_d(@dst, X, d);
    
function X = dstn_(X)
    X = make_n(@dst, X);
    
function X = idst_(X, d)
    if nargin < 2
        d = 1;
    end
    X = make_d(@idst, X, d);
    
function X = idstn_(X)
    X = make_n(@idst, X);
    
% === Helpers =============================================================

function X = make_d(f, X, d)
    if nargin < 3
        d = 1;
    end
    X = shiftdim(X, d-1);
    dim = size(X);
    X = reshape(f(reshape(X, dim(1), [])), dim);
    X = shiftdim(X, length(dim)-d+1);
    
function X = make_n(f, X)
    dim = size(X);
    n = length(dim);
    for d=1:n
        X = reshape(f(reshape(X, dim(1), [])), dim);
        X = shiftdim(X, 1);
        dim = size(X);
    end