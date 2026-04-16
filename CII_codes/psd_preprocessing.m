%% PSD PROGRESSION - SUBJECT 13AR
clear; clc; close all;
fprintf('Step 1 → Noisy, weak alpha\nStep 2 → Better, but still noisy\nStep 3 → Clean, strong alpha\n\nWe progressively improved spectral estimation by transitioning from single-channel, single-trial PSD to multi-channel occipital averaging and finally trial-averaged PSD, resulting in a stable and physiologically meaningful alpha peak.\n');
%% =========================
% LOAD DATA
%% =========================

file_name = 'C:\Users\VAIBHAV\Desktop\NSP_GrantProposal\ProjectDhyaanBK1Programs-master\data\ftDataAvgRef\013AR\EC1_ep_v8.mat'; % adjust path if needed

data = load(file_name);
data_struct = data.data;

fs = data_struct.fsample;
nTrials = length(data_struct.trial);

fprintf('Trials: %d | Fs: %.1f Hz\n', nTrials, fs);

%% =========================
% DEFINE CHANNELS
%% =========================

occipital_labels = {'O1','Oz','O2','POz','PO3','PO4','PO7','PO8','Iz'};
occipital_idx = find(ismember(data_struct.label, occipital_labels));

% fallback single channel
single_ch = find(strcmp(data_struct.label, 'Oz'));

if isempty(single_ch)
    single_ch = 1; % fallback
end

%% =========================
% STEP 1: SINGLE CHANNEL (1 TRIAL)
%% =========================

trial1 = data_struct.trial{1};

signal_single = trial1(single_ch, :);

[pxx1, f] = pwelch(signal_single, [], [], [], fs);

figure;
plot(f, 10*log10(pxx1), 'LineWidth', 2);
xlabel('Frequency (Hz)');
ylabel('Power (dB)');
title('Step 1: Single Channel (1 Trial)');
xlim([0 50]);
grid on;

%% =========================
% STEP 2: OCCIPITAL CLUSTER (1 TRIAL)
%% =========================

signal_occ = mean(trial1(occipital_idx, :), 1);

[pxx2, ~] = pwelch(signal_occ, [], [], [], fs);

figure;
plot(f, 10*log10(pxx2), 'LineWidth', 2);
xlabel('Frequency (Hz)');
ylabel('Power (dB)');
title('Step 2: Occipital Cluster (1 Trial)');
xlim([0 50]);
grid on;

%% =========================
% STEP 3: OCCIPITAL (ALL TRIALS AVG)
%% =========================

psd_all = zeros(nTrials, length(f));

for t = 1:nTrials
    trial = data_struct.trial{t};
    
    signal = mean(trial(occipital_idx,:),1);
    
    [pxx, ~] = pwelch(signal, [], [], [], fs);
    
    psd_all(t,:) = pxx;
end

psd_mean = mean(psd_all, 1);

figure;
plot(f, 10*log10(psd_mean), 'LineWidth', 2);
xlabel('Frequency (Hz)');
ylabel('Power (dB)');
title('Step 3: Occipital (All Trials Averaged)');
xlim([0 50]);
grid on;

%% =========================
% OPTIONAL: ALL IN ONE FIGURE
%% =========================

figure; hold on;

plot(f, 10*log10(pxx1), '--', 'LineWidth', 1.5);
plot(f, 10*log10(pxx2), '-', 'LineWidth', 1.5);
plot(f, 10*log10(psd_mean), 'LineWidth', 2);

legend({'Single Channel','Occipital (1 Trial)','Occipital (Avg Trials)'});
xlabel('Frequency (Hz)');
ylabel('Power (dB)');
title('PSD Progression');
xlim([0 50]);
grid on;