%% MULTI-CONDITION PSD COMPARISON - 13AR
clear; clc; close all;

%% =========================
% FILE PATHS (EDIT IF NEEDED)
%% =========================

base_path = 'C:\Users\VAIBHAV\Desktop\NSP_GrantProposal\ProjectDhyaanBK1Programs-master\data\ftDataAvgRef\013AR\';

files = {
    'EC1_ep_v8.mat'
    'EO1_ep_v8.mat'
    'G1_ep_v8.mat'
};

labels = {'Eyes Closed','Eyes Open','Grating'};

colors = lines(3);

%% =========================
% OCCIPITAL CHANNELS
%% =========================

occipital_labels = {'O1','Oz','O2','POz','PO3','PO4','PO7','PO8','Iz'};

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
    
    % Find occipital indices
    occipital_idx = find(ismember(data_struct.label, occipital_labels));
    
    psd_all = [];
    
    for t = 1:nTrials
        
        trial = data_struct.trial{t};
        
        signal = mean(trial(occipital_idx,:),1);
        
        [pxx, f] = pwelch(signal, [], [], [], fs);
        
        psd_all(t,:) = pxx;
    end
    
    % Average PSD
    psd_mean = mean(psd_all,1);
    
    % Plot
    plot(f, 10*log10(psd_mean), 'LineWidth', 2, 'Color', colors(c,:));
end

%% =========================
% PLOT SETTINGS
%% =========================

xlabel('Frequency (Hz)');
ylabel('Power (dB)');
title('PSD Comparison Across Conditions (Occipital)');
legend(labels);
xlim([0 50]);
grid on;