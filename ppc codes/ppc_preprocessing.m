%% PPC vs FREQUENCY PROGRESSION - SUBJECT 13AR
clear; clc; close all;

fprintf(['Step 1 → Single channel pair (1 trial)\n' ...
         'Step 2 → Region-level signals (1 trial)\n' ...
         'Step 3 → Region-level signals (all trials averaged)\n\n']);

%% =========================
% LOAD DATA
%% =========================

file_name = 'D:\ftavgdata\013AR\EC1_ep_v8.mat';

data = load(file_name);
data_struct = data.data;

fs = data_struct.fsample;
nTrials = length(data_struct.trial);

%% =========================
% CHANNELS
%% =========================

occipital_labels = {'O1','Oz','O2','POz','PO3','PO4','PO7','PO8','Iz'};
frontal_labels   = {'Fz','F1','F2','F3','F4'};

occ_idx = find(ismember(data_struct.label, occipital_labels));
fro_idx = find(ismember(data_struct.label, frontal_labels));

ch_occ = find(strcmp(data_struct.label,'Oz'));
ch_fro = find(strcmp(data_struct.label,'Fz'));

if isempty(ch_occ), ch_occ = occ_idx(1); end
if isempty(ch_fro), ch_fro = fro_idx(1); end

%% =========================
% FREQUENCY RANGE
%% =========================

freqs = 2:1:40;

ppc_step1 = zeros(size(freqs));
ppc_step2 = zeros(size(freqs));
ppc_step3 = zeros(size(freqs));

%% =========================
% STEP 1: SINGLE CHANNEL (1 TRIAL)
%% =========================

trial1 = data_struct.trial{1};

sig1 = trial1(ch_fro,:);
sig2 = trial1(ch_occ,:);

for i = 1:length(freqs)
    
    band = [freqs(i)-1 freqs(i)+1];
    
    ppc_step1(i) = compute_ppc_band(sig1, sig2, fs, band);
end

%% =========================
% STEP 2: REGION (1 TRIAL)
%% =========================

sig1 = mean(trial1(fro_idx,:),1);
sig2 = mean(trial1(occ_idx,:),1);

for i = 1:length(freqs)
    
    band = [freqs(i)-1 freqs(i)+1];
    
    ppc_step2(i) = compute_ppc_band(sig1, sig2, fs, band);
end

%% =========================
% STEP 3: REGION (ALL TRIALS)
%% =========================

ppc_trials = zeros(nTrials, length(freqs));

for t = 1:nTrials
    
    trial = data_struct.trial{t};
    
    sig1 = mean(trial(fro_idx,:),1);
    sig2 = mean(trial(occ_idx,:),1);
    
    for i = 1:length(freqs)
        
        band = [freqs(i)-1 freqs(i)+1];
        
        ppc_trials(t,i) = compute_ppc_band(sig1, sig2, fs, band);
    end
end

ppc_step3 = mean(ppc_trials,1);

%% =========================
% PLOT
%% =========================

figure; hold on;

plot(freqs, ppc_step1, '--', 'LineWidth',1.5);
plot(freqs, ppc_step2, '-',  'LineWidth',1.5);
plot(freqs, ppc_step3, 'LineWidth',2);

xlabel('Frequency (Hz)');
ylabel('PPC');
title('PPC vs Frequency Progression');
legend({'Single Channel','Region (1 Trial)','Region (All Trials)'});
grid on;
xlim([2 40]);

%% =========================
% FUNCTIONS
%% =========================

function ppc = compute_ppc_band(sig1, sig2, fs, band)

    % High-pass
    [b,a] = butter(3, 1/(fs/2), 'high');
    sig1 = filtfilt(b,a,sig1);
    sig2 = filtfilt(b,a,sig2);
    
    % Bandpass
    sig1 = bandpass_filt(sig1, fs, band);
    sig2 = bandpass_filt(sig2, fs, band);
    
    % Phase
    ph1 = angle(hilbert(sig1));
    ph2 = angle(hilbert(sig2));
    
    dphi = ph1 - ph2;
    
    % PPC
    N = length(dphi);
    z = exp(1i*dphi);
    ppc = (abs(sum(z))^2 - N) / (N*(N-1));
    ppc = real(ppc);
end

function y = bandpass_filt(x, fs, band)
    [b,a] = butter(3, band/(fs/2), 'bandpass');
    y = filtfilt(b,a,x);
end