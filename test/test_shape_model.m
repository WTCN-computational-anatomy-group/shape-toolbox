%function test_shape_model

N = inf;
input_dir = '/Users/balbasty/Desktop/IXI_2D/input';
files     = spm_select('List', input_dir, '\.img$');
if ~isfinite(N)
    N         = size(files, 1);
end
fnames    = cell(1, N);
for n=1:N
    fnames{n} = fullfile(input_dir, deblank(files(n,:)));
end
    
%%

N       = 10;
input   = struct;
input.f = fnames(1:N);

opt = struct;
opt.dir.model  = '/Users/balbasty/Desktop/IXI_2D/output_shape';
opt.dir.dat    = opt.dir.model;
opt.model      = struct('name', 'categorical');
opt.pg.K       = 5;
opt.lb.exact   = false;

% cluster
% opt.dist.server.ip      = '';
% opt.dist.server.login   = 'ybalba';
% opt.dist.server.folder  = '/data/ybalba/distribute';
% opt.dist.client.folder  = '/Users/balbasty/Desktop/FIL/data/distribute';
% opt.dist.matlab.bin     = '/share/apps/matlab';
% opt.dist.matlab.add     = {'/data/ybalba/matlab/shape-toolbox' ...
%                            '/data/ybalba/matlab/utility-functions'};
% opt.dist.spm.path       = '/data/ybalba/matlab/spm-trunk';
% opt.dist.spm.toolboxes  = {'Shoot'};
% opt.dist.translate      = {'/Users/balbasty/Desktop/FIL/data' '/data/ybalba'};
% opt.dist.restrict       = 'file_array';


%%
[model, dat, opt] = shape_model(input, opt);

