%% ALPHA POWER COMPARISON - SUBJECT 13AR
clear; clc; close all;

%% =========================
% PATH SETTINGS (EDIT THIS)
%% =========================

base_path = 'C:\Users\VAIBHAV\Desktop\NSP_GrantProposal\ProjectDhyaanBK1Programs-master\data\ftDataAvgRef\013AR\';

files = {
    'EC1_ep_v8.mat'   % Eyes Closed
    'EO1_ep_v8.mat'   % Eyes Open
    'G1_ep_v8.mat'   % Grating
};

labels_plot = {'EC','EO','GR'};

%% =========================
% DEFINE OCCIPITAL CHANNELS
%% =========================

occipital_labels = {'O1','Oz','O2','POz','PO3','PO4','PO7','PO8','Iz'};

%% =========================
% DEFINE ALPHA BAND
%% =========================

alpha_band = [8 12];

alpha_power = [];
condition_labels = [];

%% =========================
% MAIN LOOP
%% =========================

for c = 1:length(files)
    
    file_name = fullfile(base_path, files{c});
    
    data = load(file_name);
    data_struct = data.data;
    
    fs = data_struct.fsample;
    nTrials = length(data_struct.trial);
    
    % Find occipital channels
    occipital_idx = find(ismember(data_struct.label, occipital_labels));
    
    fprintf('Processing %s | Trials: %d\n', files{c}, nTrials);
    
    for t = 1:nTrials
        
        trial = data_struct.trial{t};
        
        % Average occipital channels
        signal = mean(trial(occipital_idx,:),1);
        
        % Compute PSD
        [pxx, f] = pwelch(signal, [], [], [], fs);
        
        % Extract alpha band power
        idx = f >= alpha_band(1) & f <= alpha_band(2);
        alpha_val = mean(pxx(idx));
        
        % Store
        alpha_power = [alpha_power; alpha_val]; %#ok<AGROW>
        condition_labels = [condition_labels; c]; %#ok<AGROW>
    end
end

%% =========================
% LOG TRANSFORM
%% =========================

alpha_power = log(alpha_power);

%% =========================
% BOXPLOT
%% =========================

figure;
boxplot(alpha_power, condition_labels);
xticklabels(labels_plot);
ylabel('Log Alpha Power');
title('Alpha Power Across Conditions (Occipital)');
grid on;

%% =========================
% PRINT MEAN VALUES
%% =========================

for c = 1:length(files)
    mean_val = mean(alpha_power(condition_labels == c));
    fprintf('Mean Alpha (%s): %.3f\n', labels_plot{c}, mean_val);
end
%% =========================
% STATISTICS
%% =========================

alpha_EC = alpha_power(condition_labels == 1);
alpha_EO = alpha_power(condition_labels == 2);
alpha_GR = alpha_power(condition_labels == 3);

% T-tests
[~, p_EC_EO] = ttest2(alpha_EC, alpha_EO);
[~, p_EC_GR] = ttest2(alpha_EC, alpha_GR);
[~, p_EO_GR] = ttest2(alpha_EO, alpha_GR);

fprintf('\nP-values:\n');
fprintf('EC vs EO: %.5f\n', p_EC_EO);
fprintf('EC vs GR: %.5f\n', p_EC_GR);
fprintf('EO vs GR: %.5f\n', p_EO_GR);