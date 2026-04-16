%% TRUE CII (ALPHA ONLY - FIXED)
clear; clc; close all;

%% PATH
base_path = 'C:\Users\VAIBHAV\Desktop\NSP_GrantProposal\ProjectDhyaanBK1Programs-master\data\ftDataAvgRef\013AR\';

files = {
    'EC1_ep_v8.mat'
    'EO1_ep_v8.mat'
    'G1_ep_v8.mat'
};

labels_plot = {'EC','EO','GR'};

%% CHANNELS
occipital_labels = {'O1','Oz','O2','POz','PO3','PO4','PO7','PO8','Iz'};
alpha_band = [8 12];

features = cell(3,1);

%% FEATURE EXTRACTION
for c = 1:3
    
    data = load(fullfile(base_path, files{c}));
    data_struct = data.data;
    
    fs = data_struct.fsample;
    occipital_idx = find(ismember(data_struct.label, occipital_labels));
    
    feats = [];
    
    for t = 1:length(data_struct.trial)
        
        trial = data_struct.trial{t};
        signal = mean(trial(occipital_idx,:),1);
        
        [pxx, f] = pwelch(signal, [], [], [], fs);
        
        idx = f >= alpha_band(1) & f <= alpha_band(2);
        feats = [feats; mean(pxx(idx))]; %#ok<AGROW>
    end
    
    feats = log(feats + 1e-8);
    features{c} = feats;
end

%% =========================
% SIMILARITY MATRIX (FIXED)
%% =========================

sim_matrix = zeros(3);

for i = 1:3
    for j = 1:3
        
        Xi = features{i};
        Xj = features{j};
        
        sims = [];
        
        for a = 1:length(Xi)
            for b = 1:length(Xj)
                
                % similarity via inverse distance
                d = abs(Xi(a) - Xj(b));
                sims = [sims; exp(-d)]; %#ok<AGROW>
            end
        end
        
        sim_matrix(i,j) = mean(sims);
    end
end

%% DISPLAY
disp('=== FIXED ALPHA SIMILARITY ===');
disp(sim_matrix);

figure;
imagesc(sim_matrix);
colorbar;
title('Alpha Similarity (Fixed)');
xticks(1:3); yticks(1:3);
xticklabels(labels_plot);
yticklabels(labels_plot);

%% CII
CII_alpha = mean(sim_matrix(~eye(3)));

fprintf('\nCII (Alpha-fixed): %.3f\n', CII_alpha);