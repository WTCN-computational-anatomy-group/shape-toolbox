%function test_pgva_model_2

N = inf;
input_dir = '/Users/balbasty/Desktop/model/input';
files     = spm_select('List', input_dir, '\.img$');
if ~isfinite(N)
    N         = size(files, 1);
end
fnames    = cell(1, N);
for n=1:N
    fnames{n} = fullfile(input_dir, deblank(files(n,:)));
end
    
%%

N       = 20;
input   = struct;
input.f = fnames(1:N);

opt = struct;
opt.dir.model = '/Users/balbasty/Desktop/model/outputvel2';
opt.dir.dat   = opt.dir.model;
opt.model     = struct('name', 'categorical');
opt.pg.K      = 9;
opt.split.par = false;

%%
[model, dat, opt] = pgva_model(input, opt);

