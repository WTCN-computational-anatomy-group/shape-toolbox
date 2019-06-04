%% Set path
% We need SPM and the Shoot toolbox to be on the path:
% https://www.fil.ion.ucl.ac.uk/spm/
% Plus these two additional toolboxes:
% https://github.com/WTCN-computational-anatomy-group/auxiliary-functions
% https://github.com/WTCN-computational-anatomy-group/distributed-computing
%
% It is advised to have them all at the same level:
% ../shape-toolbox
% ../auxiliary-functions
% ../distributed-computing
setpath;

%% Specify input/output folders

input_folder = './input/';
output_folder = './output/';

%% Download demo dataset

download_demo_set = true;
if download_demo_set
    ixi_2d_dropbox = 'https://dl.dropboxusercontent.com/s/1gf5qsbf30hd4pg/IXI_2D.zip';
    websave(fullfile(input_folder,'IXI_2D.zip'),ixi_2d_dropbox);
    unzip(fullfile(input_folder,'IXI_2D.zip'),input_folder);
end

%% Read input files

N = inf;
files     = spm_select('List', input_folder, '\.img$');
if ~isfinite(N)
    N         = size(files, 1);
end
fnames    = cell(1, N);
for n=1:N
    fnames{n} = fullfile(input_folder, deblank(files(n,:)));
end

%% Prepare input structure

% How may subjects should we use? (581 subjects in total in IXI_2D)
N       = 100;
input   = struct;
input.f = fnames(1:N);

% Specify options
opt = struct;
opt.dir.model  = output_folder;
opt.dir.dat    = opt.dir.model;
opt.model      = struct('name', 'categorical'); % Input files are segmentations
opt.pg.K       = 50;                            % Number of principal modes (should be less that N)
opt.pg.prm     = 0.1*[0.001 0 10 0.1 0.2];      % Regularisation of the velocity fields
opt.par.subjects.mode = 'parfor';               % 'for'/'parfor'


%% Train model
[model, dat, opt] = shape_model(input, opt);

