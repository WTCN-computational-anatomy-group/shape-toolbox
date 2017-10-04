clear s
s.Kv          = 48;
s.Ka          = 48;
s.linked      = true;
s.vx          = [1.5 1.5 1.5];
s.v_settings  = [1e-4  0.1   2 0.25 0.5]*0.01;
s.a_settings  = [1e-2  1e-0 0]*0.1;
s.mu_settings = [1e-4  1e-2 0];

s.mg_its      = [3 3];
s.int_args    = 8;
s.bs_args     = [1 1 0  1 1 1];
s.nit         = 2;
s.maxit       = 10;
s.omega       = 1.0;
s.wt          = [1.0 1.0];
s.likelihood  = 'multinomial';
s.result_name = 'IXIcrap';
s.result_dir  = '/Users/balbasty/Desktop/model/output';
%s.nu0         = 0.1;
s.ondisk      = true;
s.batchsize   = 16;

datdir = '/Users/balbasty/Devel/ucl/shape-toolbox-balbasty/test/data/ixi_2d';
rc1    = spm_select('List',datdir,'^sv_rc1.*img');
rc2    = spm_select('List',datdir,'^sv_rc2.*img');
N      = size(rc1,1);

if 1 
    clear dat
    for n=1:N
        dat(n).f = strvcat(fullfile(datdir,deblank(rc1(n,:))),...
                           fullfile(datdir,deblank(rc2(n,:))));
    end
else
    for n=1:N
        dat(n).f = fullfile(datdir,deblank(rc1(n,:)));
    end
    s.likelihood = 'binomial';
end

results_file = fullfile('/Users/balbasty/Desktop/model/output',['train_' s.result_name '.mat']);
PG1

