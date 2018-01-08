%function test_pg_from_vel

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
for n=1:N
    arrays{n} = datbase(n).v;
end
dat = struct('v', arrays);


opt.directory = '/Users/balbasty/Desktop/model/outputvel';
opt.vs  = optbase.vs;
opt.prm = optbase.prm;
opt.par = false;

%%
[model, dat] = pg_from_vel(opt, dat);

