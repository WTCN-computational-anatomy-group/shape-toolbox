%%% TEST SCRIPT

%% Read real velocity
vn = nifti('v_digit7_1_Template.nii');
vdata = reshape(vn.dat, [28 28 28 3]);
vel_dim =[28 28 28];
vel_vs = [1 1 1];
im_dim = 3;
param = [0.0001, 0.001, 0.2, 0.05, 0.2];
v = zeros([vel_dim im_dim]);
v(:) = vdata(:);

%% Compute L
L = diffeo('penalty', vel_dim, vel_vs, param);

%% Compute m
m = fmg(v, param, 2, 2);

%% FUNCTIONS
function x = fmg(b, param, vs, nlevels, nrelax)

end

function x = inv_relax(A, b, varargin)
% Inversion by Jacobi's relaxation method. The original method needs A to 
% be diagonally dominant. A stabilizing value is used if this is not the
% case.
% We cannot use Gauss-Siedel in Matlab because it needs iterating over
% matrix elements.
    opt = struct;
    opt.niter = 10;
    opt.safe = true;
%     opt.miniter = 0;
%     opt.maxiter = 0;
%     opt.stopval = 1E-4;
%     opt.stopcrit = @ssd;
    opt = parse_varargin(opt, varargin);
    
    % Split A into diagonal and non-diagonal parts
    dim = size(A);
    D = spdiags(diag(A), 0, dim(1), dim(2));
    E = A - D;
    
    % Find the stabilizing value if A is not diag dominant
    s = 0;
    if opt.safe
        diff = abs(diag(D)) - sum(abs(E), 2);
        if ~all(diff >= 0)
            s = - min(diff);
        end
        clear diff
    end
    
    % Initialize x
    x = ones(size(b));
    % Iterate
    for i=1:opt.niter
        if s > 0
            x = x + (1. / (D + s)) * (b - A * x);
        else
            x = (1. / D) * (b - E * x);
        end
    end
end