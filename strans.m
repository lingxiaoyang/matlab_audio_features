% Calculate features in batch.
WIN_SIZE = 2048;
HOP_SIZE = 256;
MIN_F0 = 70;
MAX_F0 = 1000;
WAV_DIR = '/highway/strans_621/waves';
FEAT_DIR = '/highway/strans_621/train';

FEATURES = {
  'SpectralCentroid',    % 1
  'SpectralCrest',       % 1
  'SpectralDecrease',    % 1
  'SpectralFlatness',    % 1
  'SpectralFlux',        % 1
  'SpectralKurtosis',    % 1
  'SpectralMfccs',       % 13
  %'SpectralPitchChroma', % 12
  'SpectralRolloff',     % 1
  'SpectralSkewness',    % 1
  'SpectralSlope',       % 1
  'SpectralSpread',      % 1
  'SpectralTonalPowerRatio',   % 1
  %%%%%'TimeAcfCoeff',
  %%%%%'TimeMaxAcf',
  'TimePeakEnvelope',    % 2
  'TimePredictivityRatio',   % 1
  'TimeRms',             % 1
  'TimeStd',             % 1
  %%%%%'TimeZeroCrossingRate',
  
  'PitchSpectralAcf',    % 1
  'PitchSpectralHps',    % 1
  'PitchTimeAcf',        % 1
  'PitchTimeAmdf',       % 1
  'PitchTimeAuditory',   % 1
  'PitchTimeZeroCrossings',    % 1
  
  'NoveltyFlux',         % 1
  'NoveltyLaroche',      % 1
  'NoveltyHainsworth',   % 1

  'YIN',  % de Cheveigne  % 3
  'pitchflow',  % https://github.com/dpwe/pitchflow   % 3
};


addpath('./plugins/yin');
addpath('./plugins/pitchflow');

%%%%%%%%%%%%%%%%%%%%%
stream_infos = [];
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
        this_stream_info = any(isnan(v), 2);
        disp(this_stream_info);
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
    this_stream_infos = any(isnan(all_features), 2);
    if i == 1
        stream_infos = this_stream_infos;
    else
        stream_infos = stream_infos | this_stream_infos;
    end
    all_features(isnan(all_features)) = -1e+10;  % hts convention
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
stream_infos
dlmwrite([FEAT_DIR '/msdInfo'], stream_infos);
