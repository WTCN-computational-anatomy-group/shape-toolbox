%function test_pgva

% load('/Volumes/NO NAME/ucl/model/output/pg_result.mat');
load('/Users/balbasty/Desktop/model/output/pg_result.mat');
%%
datbase   = dat;
modelbase = model;
optbase   = opt;
% [datbase, modelbase] = translatePath(datbase, modelbase, ...
%     'C:\Users\ybalba\', '/Volumes/NO NAME/ucl/', '/', '\');
clear dat model opt

%%
opt = struct;

N = numel(datbase);
N = 50;
arrays = cell(N,1);
input = struct;
input.v = {};
for n=1:N
    input.v{n} = datbase(n).v.fname;
end

opt = struct;
opt.dir.model = '/Users/balbasty/Desktop/model/outputvel';
opt.dir.dat   = opt.dir.model;
opt.tpl.vs    = optbase.vs;
opt.pg.K      = 9;
opt.pg.prm    = optbase.prm;
opt.split.par = false;

%%
[model, dat, opt] = pgva_model(input, opt);

