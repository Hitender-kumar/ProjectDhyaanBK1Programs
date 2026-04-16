%% VARIABILITY ANALYSIS - SUBJECT 13AR
clear; clc; close all;

%% =========================
% PATH SETTINGS (EDIT THIS)
%% =========================

base_path = 'C:\Users\VAIBHAV\Desktop\NSP_GrantProposal\ProjectDhyaanBK1Programs-master\data\ftDataAvgRef\013AR\';

files = {
    'EC1_ep_v8.mat'   % Eyes Closed
    'EO1_ep_v8.mat'   % Eyes Open
    'G1_ep_v8.mat'    % Grating
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

%% =========================
% INITIALIZE STORAGE
%% =========================

alpha_power = [];
condition_labels = [];

%% =========================
% MAIN LOOP (EXTRACT FEATURES)
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
% SEPARATE CONDITIONS
%% =========================

alpha_EC = alpha_power(condition_labels == 1);
alpha_EO = alpha_power(condition_labels == 2);
alpha_GR = alpha_power(condition_labels == 3);

%% =========================
% COMPUTE MEANS
%% =========================

mean_EC = mean(alpha_EC);
mean_EO = mean(alpha_EO);
mean_GR = mean(alpha_GR);

condition_means = [mean_EC, mean_EO, mean_GR];

%% =========================
% VARIABILITY (KEY METRIC)
%% =========================

variability = var(condition_means);

%% =========================
% DISPLAY RESULTS
%% =========================

fprintf('\n===== CONDITION MEANS =====\n');
fprintf('EC: %.3f\n', mean_EC);
fprintf('EO: %.3f\n', mean_EO);
fprintf('GR: %.3f\n', mean_GR);

fprintf('\nVariance across conditions: %.6f\n', variability);

%% =========================
% VISUALIZATION
%% =========================

figure;
bar(condition_means);
set(gca, 'XTickLabel', labels_plot);
ylabel('Mean Log Alpha Power');
title('Mean Alpha Power Across Conditions');
grid on;

%% =========================
% INTERPRETATION GUIDE (DISPLAY)
%% =========================

fprintf('\n===== INTERPRETATION =====\n');

if variability < 0.05
    fprintf('Low variability → Stable / invariant brain state\n');
elseif variability < 0.2
    fprintf('Moderate variability → Partial context dependence\n');
else
    fprintf('High variability → Strong context dependence\n');
end