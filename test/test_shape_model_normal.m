%function test_shape_model_intensity

N = inf;
input_dir = '/Users/balbasty/Desktop/OASIS-LONG-yael-ra-cr-rn-bf-ss-ni/input';
files     = spm_select('List', input_dir, '\.nii$');
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
opt.dir.model  = '/Users/balbasty/Desktop/OASIS-LONG-yael-ra-cr-rn-bf-ss-ni/output_shape';
opt.dir.dat    = opt.dir.model;
opt.model      = struct('name', 'normal', 'sigma2', 100);
opt.pg.prm     = [0.001 0 10 0.1 0.2]*1000;
opt.tpl.prm    = [0.0001 0.1 0.5]*0.01;
opt.pg.K       = 5;
opt.lb.exact   = false;
opt.optimise.q = false;

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

