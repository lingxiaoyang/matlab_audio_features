function [F0, score] = mbsc(y, fs, P)
% Multi-band summary correlogram (MBSC) pitch detector for speech.
% This algorithm is described in the paper:
% L. N. Tan, A. Alwan, "Multi-band summary correlogram-based pitch detection 
% for noisy speech", Speech Communication, in press.
% Please cite the above reference, and this MBSC code package (with the 
% SPAPL shareware URL) in your work if you make use this code.
% 
% Inputs: y - input samples
%        fs - sampling rate in Hz
%        -----vthres (optional) - threshold for voiced/unvoiced detection
%        -----                    default value = 0.375
%        P  - config struct
%
% Outputs: F0 - Estimated F0 value, frame shift = 10ms. 
%               If F0=0, it denotes an unvoiced frame. 
%          score - peak amplitude of the MBSC which denotes the frame's
%                  degree of voicing.
%
% Written by Lee Ngee Tan
% University of California (Los Angeles)
% Date: Apr 25, 2013
% Email: ngee@seas.ucla.edu
%==========================================================================

vThres = P.vThres;
minF0 = P.minf0;
maxF0 = P.maxf0;
timeshift = P.hop_size / fs;

% Resample data to 8 kHz if fs is not 8 kHz
if(fs~=8000)
    ytmp = resample(y,8000,fs);
    y = ytmp;
    fs = 8000;
    clear ytmp;
end

timeshift = 0.01;
winshift = timeshift*fs;
if (maxF0>500)
    disp(['Maximum F0 should not exceed 500 Hz because multiple harmonics',...
        'will not be captured by the 1 kHz subband bandwidth']);
    return;
end
% A maximum of 500 Hz is allowed for the current design.
minlag = floor(fs/maxF0);
maxlag = ceil(fs/minF0);
n_period = 4.0;
n_sb = 4;

y_len = length(y);
n_fr = floor(y_len/winshift);
   
F0cand = zeros(n_fr,10);
cand_scores = zeros(n_fr,10);  
F0 = zeros(n_fr,1);
score = zeros(n_fr,1);  

lag = minlag-1:maxlag+1;
doublag = (maxlag+2):2:2*maxlag;
lagRange = [lag, doublag];
n_ch = length(lagRange);
quantWinlen = 80:40:640;
winlen = zeros(1,n_ch);
% Generate quantized window lengths
for k = 1:n_ch  
    [minval,minidx] = min(abs(n_period*lagRange(k)-quantWinlen));
    winlen(k) = quantWinlen(minidx);
end
halfwinlen = floor(winlen/2);
Rnfft = zeros(n_ch,1);
nfft = 8192;
for ch = 1:n_ch
    Rnfft(ch) = 2^(nextpow2(winlen(ch))+1);
end

n_stream = n_sb+1;
sb_ch_snr = zeros(n_stream,n_ch);
smooth_sb_ch_snr = zeros(n_stream,n_ch);
minSNRpk = 1.0;
sR_smooth = zeros(1,maxlag+1);
lagwt_line = polyfit([1,maxlag],[1,0.3],1);
sb_y = zeros(y_len,n_stream);
halfFs = fs/2;
n_firpts = 32;
sb_fir = zeros(n_sb,n_firpts+1);
sb_fir(1,:) = fir1(n_firpts, 1000/halfFs);
sb_fir(2,:) = fir1(n_firpts,[800, 1800]/halfFs);
sb_fir(3,:) = fir1(n_firpts,[1600, 2600]/halfFs);
sb_fir(4,:) = fir1(n_firpts,[2400, 3400]/halfFs);

% Generate signal and noise-capturing comb filters
[combFiltmat, noise_combFiltmat, expdecaycell, n_posfft]...
                            = gen_LPcombFBKs(lagRange, nfft, fs);

for sb = 1:n_sb
    tmp = filter(sb_fir(sb,:),1,[y;zeros(n_firpts/2,1)]);
    if(sb==1)
        sb_y(:,sb) = tmp(n_firpts/2+1:end);
    end
    sb_y(:,sb+1) = abs(hilbert(tmp(n_firpts/2+1:end))).^2;   
end

lim1 = find(lagRange==minlag-1);
lim2 = find(lagRange==maxlag+1);

combfilt_spec = zeros(n_posfft,n_stream,n_ch);
sb_timerel_ctr = ones(1,n_stream);
        
for fr = 1:n_fr % for each time frame
    % Perform framing and mean normalization
    for ch = 1:n_ch % for each channel 
        if(ch>1)
            if(winlen(ch)~=winlen(ch-1))
                j_start = fr*winshift+1-halfwinlen(ch);
                j_end = j_start+winlen(ch);
                ybuf4fft = zeros(winlen(ch)+1,n_stream);

                for sb = 1:n_stream
                    if(sb==1)
                        if(j_start<1)
                            sb_z_fr = [zeros(abs(j_start)+1,1);sb_y(1:j_end,sb)];
                        elseif(j_end>y_len)
                            sb_z_fr = [sb_y(j_start:y_len,sb);zeros(j_end-y_len,1)];
                        else
                            sb_z_fr = sb_y(j_start:j_end,sb);
                        end
                    else                        
                        if(j_start<1)
                            sb_z_fr = [zeros(abs(j_start)+1,1);remove_mean(sb_y(1:j_end,sb))];
                        elseif(j_end>y_len)
                            sb_z_fr = [remove_mean(sb_y(j_start:y_len,sb));zeros(j_end-y_len,1)];
                        else
                            sb_z_fr = remove_mean(sb_y(j_start:j_end,sb));
                        end
                    end
                    ybuf4fft(1:winlen(ch)+1,sb) = sb_z_fr;
                end
                fr_spec = fft(ybuf4fft,nfft);
         
            end
            
        else
            j_start = fr*winshift+1-halfwinlen(ch);
            j_end = j_start+winlen(ch);
            ybuf4fft = zeros(winlen(ch)+1,n_stream);

            for sb = 1:n_stream
                if(sb==1)
                    if(j_start<1)
                        sb_z_fr = [zeros(abs(j_start)+1,1);sb_y(1:j_end,sb)];
                    elseif(j_end>y_len)
                        sb_z_fr = [sb_y(j_start:y_len,sb);zeros(j_end-y_len,1)];
                    else
                        sb_z_fr = sb_y(j_start:j_end,sb);
                    end
                else                        
                    if(j_start<1)
                        sb_z_fr = [zeros(abs(j_start)+1,1);remove_mean(sb_y(1:j_end,sb))];
                    elseif(j_end>y_len)
                        sb_z_fr = [remove_mean(sb_y(j_start:y_len,sb));zeros(j_end-y_len,1)];
                    else
                        sb_z_fr = remove_mean(sb_y(j_start:j_end,sb));
                    end
                end
                ybuf4fft(1:winlen(ch)+1,sb) = sb_z_fr;
            end
            fr_spec = fft(ybuf4fft,nfft);
            
        end
                
        for sb = 1:n_stream  % for each subband stream
            % Perform comb filtering
            if(sb==1)
                combfilt_spec(:,sb,ch) = fr_spec(1:n_posfft,sb).*combFiltmat(:,ch).*expdecaycell(:,ch);
                noise_combfilt_spec = fr_spec(1:n_posfft,sb).*noise_combFiltmat(:,ch).*expdecaycell(:,ch);
            else
                combfilt_spec(:,sb,ch) = fr_spec(1:n_posfft,sb).*combFiltmat(:,ch);
                noise_combfilt_spec = fr_spec(1:n_posfft,sb).*noise_combFiltmat(:,ch);
            end

            % Calculate HSR (pseudo-SNR)
            sb_ch_sigpow = combfilt_spec(:,sb,ch)'*combfilt_spec(:,sb,ch);
            sb_ch_noisepow = noise_combfilt_spec'*noise_combfilt_spec;
            sb_ch_snr(sb,ch) = sb_ch_sigpow/sb_ch_noisepow;

        end
    end
    
    sb_normR = zeros(n_ch,maxlag+1,n_stream);
    smooth_snrWt_normR = zeros(n_stream,maxlag+1);
            
    sb_snr_wt = zeros(n_stream,n_ch);
    smooth_sb_snr_wt = zeros(n_stream,n_ch);
    smooth_max_chan_snr = zeros(1,n_stream);
        
    for sb = 1:n_stream
        if(fr==1) % 1st update
            smooth_sb_ch_snr(sb,:) = sb_ch_snr(sb,:);
        else % Time-smoothing of HSR
            smooth_sb_ch_snr(sb,:) = 0.5*(sb_ch_snr(sb,:) + smooth_sb_ch_snr(sb,:));
        end


        % Perform channel selection
        % Stage 1 channel selection
        [snrpk,snrpkloc] = findpeaks_fast(smooth_sb_ch_snr(sb,lim1:end),minSNRpk,'nosort');           
        snrpkloc = snrpkloc + lim1 - 1;
        % Stage 2 channel selection
        if(~isempty(find(snrpkloc<lim2, 1)))
            snrpklag = lagRange(snrpkloc);
            n_pks = length(snrpklag);
            validPks = false(n_pks,1);
            for pk = 1:n_pks
                if(snrpklag(pk)<=maxlag)
                    doublagloc = find(abs(0.5*snrpklag/snrpklag(pk)-1)<0.2);
                    if(~isempty(doublagloc))
                        if(length(doublagloc)>1)
                            [min_val,min_idx] = min(abs(0.5*snrpklag(doublagloc)/snrpklag(pk)-1));
                            validPks([pk,doublagloc(min_idx)])=1;
                        else
                            validPks([pk,doublagloc])=1;
                        end
                    end
                else
                    break;
                end
            end
            valid_snrpklocs = snrpkloc(validPks);
            if(sum(validPks))
                smooth_sb_snr_wt(sb,valid_snrpklocs) = smooth_sb_ch_snr(sb,valid_snrpklocs)-1;
            end
        end    

        pkChIdx = find(smooth_sb_snr_wt(sb,:)>0);
        % Stage 3 channel selection
        for pkch = 1:length(pkChIdx)
            ch = pkChIdx(pkch);

            % Compute energy-normalized autocorrelation for selected channels
            n_Rlag = min(winlen(ch)-1,maxlag);

            tmp = ifft([combfilt_spec(:,sb,ch);zeros(nfft-2*n_posfft+1,1);conj(combfilt_spec(end:-1:2,sb,ch))],nfft);
            sb_comb_y = tmp(1:winlen(ch));
            tmp_sb_normR = normAC(sb_comb_y, min(winlen(ch)-1,lagRange(end)));
            sb_normR(ch,1:n_Rlag+1,sb) = tmp_sb_normR(1:n_Rlag+1);
            [maxpk,pklag] = findMaxPk(tmp_sb_normR);
            if(~isempty(pklag))
                if(abs(1-pklag/lagRange(ch))>0.2 && abs(1-2*pklag/lagRange(ch))>0.2)
                    sb_snr_wt(sb,ch) = 0;
                    smooth_sb_snr_wt(sb,ch) = 0;
                end
            else
                sb_snr_wt(sb,ch) = 0;
                smooth_sb_snr_wt(sb,ch) = 0;
            end
        end
        % Compute subband summary correlogram 
        sum_smooth_sb_snr_wt = sum(smooth_sb_snr_wt(sb,:));
        if(sum_smooth_sb_snr_wt>0)
            maxpk = max(smooth_sb_snr_wt(sb,:));
            smooth_max_chan_snr(sb) = maxpk;
            smooth_sb_snr_wt(sb,:) = smooth_sb_snr_wt(sb,:)/sum_smooth_sb_snr_wt;
            
        end
        smooth_snrWt_normR(sb,:) = smooth_sb_snr_wt(sb,:)*sb_normR(:,:,sb); 
    end
    
    % Calculate between-subband reliability
    smooth_ACmaxpk_lag = zeros(1,n_stream);
    smooth_sb_rel_ctr = ones(1,n_stream);
    for sb = 1:n_stream
        if(smooth_snrWt_normR(sb,1))
            [maxpk,pkIdx] = findMaxPk(smooth_snrWt_normR(sb,:));
            if(~isempty(pkIdx))
                smooth_ACmaxpk_lag(sb) = pkIdx;
            end
        end
    end
    for sb = 1:n_stream
        lag_ratio = smooth_ACmaxpk_lag/smooth_ACmaxpk_lag(sb);
        smooth_sb_rel_ctr(sb) = sum(abs(1-lag_ratio)<0.1);
        if(fr>1)
            lag_ratio = smooth_ACmaxpk_lag(sb)/prev_ACmaxpk_lag(sb);
            if(abs(1-lag_ratio)<0.1)
                sb_timerel_ctr(sb) = sb_timerel_ctr(sb)+1;
            else
                sb_timerel_ctr(sb) = 1;
            end
        else
            sb_timerel_ctr(sb) = 1;
        end
    end
    
    prev_ACmaxpk_lag = smooth_ACmaxpk_lag;
    
    
    % Compute subband reliability & multi-band summary correlogram 
    if(sum(smooth_max_chan_snr))
        smooth_max_chan_snr2 = smooth_max_chan_snr.*smooth_sb_rel_ctr;
        smooth_max_chansnr_wt = smooth_max_chan_snr2/sum(smooth_max_chan_snr2);
        sR_smooth_sb_snr = smooth_max_chansnr_wt*smooth_snrWt_normR;        
    else
        sR_smooth_sb_snr = zeros(1,maxlag+1);
    end
    
    % Time-smooth MBSC
    if(sR_smooth(1)==0) % 1st update of sR_smooth 
        sR_smooth = sR_smooth_sb_snr;
    else
        sR_smooth = 0.5*(sR_smooth + sR_smooth_sb_snr);
    end
    
    % Peak extraction and lag-weighting
    [pks, pklocs] = findpeaks_fast(sR_smooth(minlag-1:maxlag+1),0,'descend');
    if(~isempty(pks))
        pklocs = pklocs + minlag - 2;
        len = length(pks);
        n_cands = min(len,10);
        lag_cands = zeros(n_cands,1);
        curr_F0cand = zeros(n_cands,1);
        curr_scores = zeros(n_cands,1);
        for m = 1:n_cands
            fit_x = pklocs(m)+[-1,0,1];
            p_coeff = polyfit(fit_x,sR_smooth(fit_x),2);
            x_turnPt = -0.5*p_coeff(2)/p_coeff(1);
            y_turnPt = polyval(p_coeff,x_turnPt); 
            lag_cands(m) = x_turnPt;
            curr_F0cand(m) = fs./x_turnPt;                
            curr_scores(m) = y_turnPt;
        end
        
        F0cand(fr,1:n_cands) = curr_F0cand';
        cand_scores(fr,1:n_cands) = curr_scores';
    end
         
end

% Perform V/UV detection and F0 candidate selection
vuv = cand_scores(:,1)>=vThres;
vuv = medfilt1(vuv,5);
for fr = 1:n_fr
    cands = cand_scores(fr,:)>0;
    if(sum(cands))
        lags = (fs./F0cand(fr,cands));
        lagwt = polyval(lagwt_line,lags);
        wt_scores = cand_scores(fr,cands).*lagwt;
        [max_wt_score, max_idx] = max(wt_scores);
        F0(fr) = F0cand(fr,max_idx(1));
        score(fr) = cand_scores(fr,max_idx(1));
    end
end
F0 = F0.*vuv;
    
        
% Generate signal and noise-capturing comb filterbanks (FBK)
function [combFiltmat,noise_combFiltmat, expdecaycell, n_posfft] = gen_LPcombFBKs(all_lags, nfft, Fs)

lp_freq = 1000;
fq = Fs./all_lags;
n_fq = length(fq); % Total # channels in each FBK
halfFs = Fs/2;
freq = (-halfFs+Fs/nfft:Fs/nfft:halfFs)';
freq0 = freq(nfft/2:end);
lowfcut_startIdxH = find(freq0>lp_freq, 1);
n_posfft = lowfcut_startIdxH;
combFiltmat = zeros(n_posfft,n_fq); % comb_low
noise_combFiltmat = zeros(n_posfft,n_fq); % noise_comb_low
expdecaycell = zeros(n_posfft,n_fq); % harmonic decreasing function
for ch = 1:n_fq      
    filt = 0.5*(1+cos(2*pi/fq(ch)*freq0));
    [lowfcutL,lowfcut_idxL] = min(filt(freq0/fq(ch)>=0 & freq0/fq(ch)<1)); 
    combFiltmat((lowfcut_idxL:lowfcut_startIdxH),ch) = filt(lowfcut_idxL:lowfcut_startIdxH);
    noise_combFiltmat((lowfcut_idxL:lowfcut_startIdxH),ch) = (1-filt(lowfcut_idxL:lowfcut_startIdxH));
    expdecaycell(:,ch) = 1+[1;min(1,(fq(ch)./freq0(2:n_posfft)))];
end

function y_mean0 = remove_mean(y)
y_mean0 = y - mean(y);

function c = normAC(y, maxlag)
u1 = mean(y);
y = y - u1;
L = length(y);
Rnfft = 2^(nextpow2(L)+1);
tmpf = fft(y,Rnfft);
tmpR = ifft(tmpf.*conj(tmpf),Rnfft);
c = zeros(1,maxlag+1);
c(1:min(L,maxlag+1)) = tmpR(2:min(L+1,maxlag+2))/tmpR(1);

function [pkMag, pkIdx] = findMaxPk(x)
x = x(:);
xL = x(2:end-1)>x(1:end-2);
xR = x(2:end-1)>x(3:end);
pklocs = find([0; xL&xR]==1);
pks = x(pklocs);
[pkMag, idx] = max(pks);
pkIdx = pklocs(idx);

function [pks, pklocs] = findpeaks_fast(x,minPkMag,sortstr)
x = x(:);
xL = x(2:end-1)>x(1:end-2);
xR = x(2:end-1)>x(3:end);
pklocs = find([0; xL&xR]==1);
pks = x(pklocs);
retain_pkIdx = pks>=minPkMag;
pks = pks(retain_pkIdx);
pklocs = pklocs(retain_pkIdx);
if(~isempty(pks))
    if(strcmp(sortstr,'descend'))
        [sortPks,sortIdx] = sort(pks,'descend');   
        pks = sortPks;
        pklocs = pklocs(sortIdx); 
    elseif(strcmp(sortstr,'ascend'))
        [sortPks,sortIdx] = sort(pks);
        pks = sortPks;
        pklocs = pklocs(sortIdx); 
    end
end
