%% MULTI-CONDITION PPC COMPARISON - 13AR
clear; clc; close all;

%% =========================
% FILE PATHS
%% =========================

base_path = 'D:\ftavgdata\013AR';

files = {
    'EC1_ep_v8.mat'
    'EO1_ep_v8.mat'
    'G1_ep_v8.mat'
};

labels = {'Eyes Closed','Eyes Open','Grating'};
colors = lines(3);

%% =========================
% CHANNELS
%% =========================

occipital_labels = {'O1','Oz','O2','POz','PO3','PO4','PO7','PO8','Iz'};
frontal_labels   = {'Fz','F1','F2','F3','F4'};

%% =========================
% FREQUENCY BANDS
%% =========================

bands = [
    4 8;    % theta
    8 12;   % alpha
    13 30   % beta
];

band_labels = {'Theta','Alpha','Beta'};

figure; hold on;

%% =========================
% LOOP OVER CONDITIONS
%% =========================

for c = 1:length(files)
    
    file_name = fullfile(base_path, files{c});
    
    data = load(file_name);
    data_struct = data.data;
    
    fs = data_struct.fsample;
    nTrials = length(data_struct.trial);
    
    % Channel indices
    occ_idx = find(ismember(data_struct.label, occipital_labels));
    fro_idx = find(ismember(data_struct.label, frontal_labels));
    
    ppc_all = zeros(nTrials, size(bands,1));
    
    for t = 1:nTrials
        
        trial = data_struct.trial{t};
        
        sig_occ = mean(trial(occ_idx,:),1);
        sig_fro = mean(trial(fro_idx,:),1);
        
        % High-pass
        [b,a] = butter(3, 1/(fs/2), 'high');
        sig_occ = filtfilt(b,a,sig_occ);
        sig_fro = filtfilt(b,a,sig_fro);
        
        for bnd = 1:size(bands,1)
            
            band = bands(bnd,:);
            
            % Bandpass
            s1 = bandpass_filt(sig_fro, fs, band);
            s2 = bandpass_filt(sig_occ, fs, band);
            
            % Phase
            ph1 = angle(hilbert(s1));
            ph2 = angle(hilbert(s2));
            
            dphi = ph1 - ph2;
            
            % PPC
            ppc_all(t, bnd) = compute_ppc(dphi);
        end
    end
    
    % Average across trials
    ppc_mean = mean(ppc_all,1);
    
    % Plot (like PSD but discrete)
    plot(1:length(bands), ppc_mean, '-o', ...
        'LineWidth',2, ...
        'Color',colors(c,:), ...
        'MarkerSize',8);
end

xticks(1:length(bands));
xticklabels(band_labels);

xlabel('Frequency Bands');
ylabel('PPC');
title('PPC Comparison Across Conditions (Frontal–Occipital)');
legend(labels);
grid on;

%% =========================
% FUNCTIONS
%% =========================

function y = bandpass_filt(x, fs, band)
    [b,a] = butter(3, band/(fs/2), 'bandpass');
    y = filtfilt(b,a,x);
end

function ppc = compute_ppc(dphi)
    dphi = dphi(:);
    N = length(dphi);
    z = exp(1i * dphi);
    ppc = (abs(sum(z))^2 - N) / (N*(N-1));
    ppc = real(ppc);
end
%% =========================
% PPC vs FREQUENCY (PSD-LIKE VISUALIZATION)
%% =========================

freqs = 2:1:40;   % frequency range
figure; hold on;

for c = 1:length(files)
    
    file_name = fullfile(base_path, files{c});
    
    data = load(file_name);
    data_struct = data.data;
    
    fs = data_struct.fsample;
    nTrials = length(data_struct.trial);
    
    % Channel indices
    occ_idx = find(ismember(data_struct.label, occipital_labels));
    fro_idx = find(ismember(data_struct.label, frontal_labels));
    
    ppc_trials = zeros(nTrials, length(freqs));
    
    for t = 1:nTrials
        
        trial = data_struct.trial{t};
        
        sig_occ = mean(trial(occ_idx,:),1);
        sig_fro = mean(trial(fro_idx,:),1);
        
        % High-pass
        [b,a] = butter(3, 1/(fs/2), 'high');
        sig_occ = filtfilt(b,a,sig_occ);
        sig_fro = filtfilt(b,a,sig_fro);
        
        for fidx = 1:length(freqs)
            
            band = [freqs(fidx)-1 freqs(fidx)+1];
            
            s1 = bandpass_filt(sig_fro, fs, band);
            s2 = bandpass_filt(sig_occ, fs, band);
            
            ph1 = angle(hilbert(s1));
            ph2 = angle(hilbert(s2));
            
            dphi = ph1 - ph2;
            
            ppc_trials(t, fidx) = compute_ppc(dphi);
        end
    end
    
    % Average across trials
    ppc_mean_curve = mean(ppc_trials,1);
    
    % Plot
    plot(freqs, ppc_mean_curve, ...
        'LineWidth',2, ...
        'Color',colors(c,:));
end

xlabel('Frequency (Hz)');
ylabel('PPC');
title('PPC vs Frequency (Frontal–Occipital)');
legend(labels);
grid on;
xlim([2 40]);