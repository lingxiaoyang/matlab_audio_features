% Calculate features in batch.
WIN_SIZE = 2048;
HOP_SIZE = 512;
MIN_F0 = 70;
MAX_F0 = 900;
WAV_DIR = '/highway/strans_621/waves';
FEAT_DIR = '/highway/strans_621/train';

FEATURES = {
  'SpectralCentroid',
  'SpectralCrest',
  'SpectralDecrease',
  'SpectralFlatness',
  'SpectralFlux',
  'SpectralKurtosis',
  %%'SpectralMfccs',
  %%%%%'SpectralPitchChroma',
  'SpectralRolloff',
  'SpectralSkewness',
  'SpectralSlope',
  'SpectralSpread',
  'SpectralTonalPowerRatio',
  %%%%%'TimeAcfCoeff',
  %%%%%'TimeMaxAcf',
  'TimePeakEnvelope',
  'TimePredictivityRatio',
  'TimeRms',
  'TimeStd',
  %%%%%'TimeZeroCrossingRate',
  
  'PitchSpectralAcf',
  'PitchSpectralHps',
  'PitchTimeAcf',
  'PitchTimeAmdf',
  'PitchTimeAuditory',
  'PitchTimeZeroCrossings',
  
  'NoveltyFlux',
  'NoveltyLaroche',
  'NoveltyHainsworth',

  'YIN',  % de Cheveigne
  'pitchflow',  % https://github.com/dpwe/pitchflow
};


addpath('./plugins/yin');
addpath('./plugins/pitchflow');

%%%%%%%%%%%%%%%%%%%%%
wavfiles = dir([WAV_DIR '/*.wav']);
for i = 1:length(wavfiles)
    wavfile = wavfiles(i);
    [path, base, ext] = fileparts(wavfile.name);
    [audiodata, fs] = wavread([WAV_DIR '/' wavfile.name]);
    
    % calculate feat ASCII
    featfilename = [FEAT_DIR '/' base '.feat_ascii'];
    disp(wavfile.name);
    if fs ~= 44100
        disp(['WARNING: fs should be 44100, instead of: ', num2str(fs), ' in ', wavfile.name]);
        continue;
    end


    all_features = [];
    for j = 1:length(FEATURES)
        feat_name = FEATURES{j};
        if strcmp(feat_name, 'YIN')
            p.minf0 = MIN_F0;
            p.maxf0 = MAX_F0;
            p.thresh = 0.1;
            p.hop = HOP_SIZE;
            p.bufsize = WIN_SIZE;
            p.sr = fs;
            r = yin(audiodata, p);
            v = [r.f0; r.ap; r.pwr];
            v(isnan(v)) = 0;
        elseif strcmp(feat_name, 'pitchflow')
            p.t_win = WIN_SIZE / fs;
            p.t_hop = HOP_SIZE / fs;
            vall = pitchflow(audiodata, fs, p);
            v = pitchflow_collapse(vall);
        elseif strcmp(feat_name(1:5), 'Pitch')
            name = feat_name(6:end);
            [f, t] = ComputePitch(name, audiodata, fs, [], WIN_SIZE, HOP_SIZE);
            v = f;   % f -- seems to be in log domain
        elseif strcmp(feat_name(1:7), 'Novelty')
            name = feat_name(8:end);
            [d, t, iPeaks] = ComputeNoveltyFunction(name, audiodata, fs, [], WIN_SIZE, HOP_SIZE);
            v = d;
        else
            [v, t] = ComputeFeature(FEATURES{j}, audiodata, fs, [], WIN_SIZE, HOP_SIZE);
        end
        
        dim = size(v, 1);
        disp(['INFO: dimension of feature ' FEATURES{j} ': ' num2str(dim) ' length: ' num2str(size(v, 2))]);

        laf = size(all_features, 2);
        lv = size(v, 2);
        if laf > lv && laf ~= 0
            all_features = [all_features(:,1:lv); v];
        elseif laf < lv && laf ~= 0
            all_features = [all_features; v(:, 1:laf)];
        else
            all_features = [all_features; v];
        end
    end
    
    % calculate delta and delta-delta
    dims = size(all_features, 1);
    delta_all_features = [zeros(dims, 1) diff(all_features, 1, 2)];
    len = size(all_features, 2);
    delta_delta_all_features = zeros(dims, len);
    delta_delta_all_features(:,1) = all_features(:,2)-all_features(:,1);
    for j = 2:(len-1)
        delta_delta_all_features(:,j) = (all_features(:,j+1)-all_features(:,j-1))/2;
    end
    delta_delta_all_features(:,len) = (all_features(:,len)-all_features(:,len-1));

    all_features = [all_features;delta_all_features;delta_delta_all_features];
    all_features(isnan(all_features)) = 0;
    dims = size(all_features, 1);
    size(all_features)
    
    all_features_reshaped = reshape(all_features, [], 1);
    dlmwrite(featfilename, all_features_reshaped);

    % x2x & htkheader
    cmd = ['/usr/local/SPTK/bin/x2x +af ' FEAT_DIR '/' base '.feat_ascii > ' FEAT_DIR '/' base '.feat_f'];
    disp(cmd);
    system(cmd);

    byte_per_frame = 4 * dims;
    cmd = ['/highway/strans_621/addhtkheader.pl ' num2str(fs) ' ' ...
        num2str(HOP_SIZE) ' '  num2str(byte_per_frame) ' 9 '...
        FEAT_DIR '/' base '.feat_f > ' ...
        FEAT_DIR '/' base '.feat'];
    disp(cmd);
    system(cmd);
    %%break;
end

