%% PITCHFLOW - direct calculation of delta-pitch
%
% Pitch is helpful in speech recogntion, but most often it is not
% the absolute pitch that is used (since this varies a lot by
% speaker), but some locally-normalized pitch contour.  In fact, it
% may be that only the pitch derivative is important.  Several
% researchers have looked at calculating local pitch derivatives
% without having to go through the laborious and error-prone
% process of first finding a pitch contour.  This code takes that
% approach.  
%
% The algorithm works by calculating a spectrogram of the speech,
% then warping the frequency axis to be logarithmic.  In this
% projection, a change in pitch corresponds to a simple translation
% of the position of all harmonics.  Thus, delta-pitch can be
% identified with simple cross-correlation between successive
% short-time spectra.

%% Example

% Load a sound file
[d,sr] = audioread(['/u/drspeech/data/swordfish/code/ehist/' ...
                    'BABEL_OP1_206_65882_20121201_174526_outLine.sph']);

% Plot its log-frequency spectrogram
subplot(311)
logfsgram(d,256,sr);
caxis([-30 30]);

% Calculate the normalized cross-correlation of adjacent log-f spectra
sxc = pitchflow(d, sr);

subplot(312)
tt = [0:(size(sxc,2)-1)]*0.010;
rr = [1:size(sxc,1)] - (size(sxc, 1)+1)/2;
imagesc(tt, rr, sxc); axis('xy');
grid

% Collapse each excerpt from the cross-corrletions into 3 features,
% the first three moments of the exponentiated NCC
dpf = pitchflow_collapse(sxc);

subplot(313)
plot(tt, dpf);

% Overplot the first moment on the xcorr, to show it tracks the
% main peak
subplot(312)
hold on; 
plot(tt, dpf(2,:), '-w'); 
hold off

% Line up and zoom in
linkaxes([subplot(311),subplot(312),subplot(313)], 'x')
axis([196 200 -5 10])


%% Bulk calculation
%
% Say you want to calculate this feature for a whole directory full
% of Babel utterances.  Here's how:

corpus = 'BABEL_OP1_102_LLP';
babelcorproot = '/u/drspeech/data/swordfish/corpora';
wavdevdir = fullfile(babelcorproot, corpus, 'conversational/dev/audio');
wavtrndir = fullfile(babelcorproot, corpus, 'conversational/training/audio');
ftrdevdir = fullfile(corpus, 'dev/dpitch');
ftrtrndir = fullfile(corpus, 'training/dpitch');
mymkdir(ftrdevdir);
mymkdir(ftrtrndir);
pitchflow_processdir(wavtrndir, ftrtrndir);
pitchflow_processdir(wavdevdir, ftrdevdir);
ftr2devdir = fullfile(corpus, 'dev/dpflow');
ftr2trndir = fullfile(corpus, 'training/dpflow');
mymkdir(ftr2devdir);
mymkdir(ftr2trndir);
pitchflow_reprocessdir(ftrdevdir, ftr2devdir);
pitchflow_reprocessdir(ftrtrndir, ftr2trndir);

%% Python Port
%
% The full pitchflow feature calculation pipeline has been ported
% to Python.  See the
% <https://github.com/dpwe/pitchflow pitchflow> package on GitHub.

%% Changelog
%

% 2014-02-11 v0.1  Cleaned up and added full "matlab publish" output
%
% 2014-01-28 v0.0  Initial release
%

%% Acknowledgment
%
% This work was supported by IARPA under the Babel program via a 
% subcontract from the ICSI-led team Swordfish
