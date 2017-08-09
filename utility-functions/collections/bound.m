function varargout = bound(id, varargin)
% FORMAT [j, (p)] = bound('circulant', i, m)
% FORMAT [j, (p)] = bound('neuman', i, m)
%
% Boundary tools.
%
% FORMAT help bound>function
% Returns the help file of the selected function.
    switch id
        case 'circulant'
            [varargout{1:nargout}] = circulant(varargin{:});
        case 'neuman'
            [varargout{1:nargout}] = neuman(varargin{:});
        otherwise
            error('No function named %s\n', string(id));
    end

% === Subfunctions ========================================================

function [i, p] = circulant(i, m)
% FORMAT [j, (p)] = bound('circulant', i, m)
% i - (infinite signed) index, i.e. virtual index in (-inf:1:inf)
% m - Length of the true domain.
% j - Correspondig index in the true domain with circulant conditions
% p - Signal period to which belongs the virtual index.
%
% Circulant boundary conditions: f(-1) = f(m-1) // f(m) = f(0).
% /!\ Here, indexing is C-style (true domain is [0, m-1]).
%     Subtract 1 to the input and add 1 to the output to use with matlab.

    p = floor(i/m);
    if i >= 0
        i = mod(i, m);
    else
        i = mod( m + mod(i,m), m);
    end
    
function [i, p] = neuman(i, m)
% FORMAT [j, (p)] = bound('neuman', i, m)
% i - (infinite signed) index, i.e. virtual index in (-inf:1:inf)
% m - Length of the true domain.
% j - Correspondig index in the true domain with Neuman conditions
% p - Signal period to which belongs the virtual index.
%
% Neuman boundary conditions: f(-1) = f(0) // f(m) = f(m-1).
% Equivalent to type II mirroring conditions (mirror2).
% /!\ Here, indexing is C-style (true domain is [0, m-1]).
%     Subtract 1 to the input and add 1 to the output to use with matlab.

    p = floor(i/m);
    m2 = m*2;
    if i >= 0
        i = mod(i, m2);
    else
        i = m2 - mod(-i-1, m2) - 1;
    end
    if i >= m
        i = m2 - i - 1;
    end
    