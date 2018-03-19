%function test_pgva_model_2

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

N       = 100;
input   = struct;
input.f = fnames(1:N);

opt = struct;
opt.dir.model  = '/Users/balbasty/Desktop/IXI_2D/outputvel2';
opt.dir.dat    = opt.dir.model;
opt.model      = struct('name', 'categorical');
opt.pg.K       = 19;
opt.split.par  = 0;

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
[model, dat, opt] = pgva_model(input, opt);

