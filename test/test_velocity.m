%% Load data
path = which('test_velocity');
if isempty(path)
    error(['The folder containing this function must be added to the ' ...
           'path so that I can locate test data.'])
end
path = fullfile(fileparts(path), 'data\train_mnist.mat');
load(path)
%% Initialize objects
v               = Velocity('Directory', 'D:\output');
v.MatchingTerm  = 'binomial';
v.RegParam      = s.v_settings;
% v.Mu            = mu2;
v.LogMu         = mu;
v.Image         = dat(1).f;
v.W             = Wv;
v.SigmaR        = 1;
v.Interpolation = [1 1 0 1 1 1];
v.Verbose       = true;
v.MaxGNIt       = 20;
% v.MaxLSIt       = 20;
v.Debug         = false;
%% Registration: full velocity optimisation
v.initZ('discard');
v.initR();
v.updateR();
%% Registration: latent coordinates optimisation
v.initR('discard');
v.initZ();
v.updateZ();
%% Registration: latent coordinates + residual field
% You can try with different v.SigmaR
% > If it is too high, the residual field takes over and we're stuck in a 
%   local optimum.
% > If it's low enough, the shape model takes over and converges better and
%   faster.
v.initZ();
v.initR();
v.update();
%% Look at uncertainty estimates
% v.uncertaintyR(); % < That's bad because only likelihood term
v.uncertaintyZ();
v.uncertaintyWZ();
% v.uncertaintyVel(); % < That's bad because sR is bad