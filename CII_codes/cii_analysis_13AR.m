%% TRUE CII USING TRIAL-LEVEL SIMILARITY
clear; clc; close all;

%% =========================
% PATH
%% =========================

base_path = 'C:\Users\VAIBHAV\Desktop\NSP_GrantProposal\ProjectDhyaanBK1Programs-master\data\ftDataAvgRef\013AR\';

files = {
    'EC1_ep_v8.mat'
    'EO1_ep_v8.mat'
    'G1_ep_v8.mat'
};

labels_plot = {'EC','EO','GR'};

%% =========================
% CHANNELS + BANDS
%% =========================

occipital_labels = {'O1','Oz','O2','POz','PO3','PO4','PO7','PO8','Iz'};

bands = [
    1 4;
    4 8;
    8 12;
    13 30
];

nBands = size(bands,1);

%% =========================
% FEATURE EXTRACTION (TRIAL LEVEL)
%% =========================

features = cell(3,1);

for c = 1:3
    
    file_name = fullfile(base_path, files{c});
    data = load(file_name);
    data_struct = data.data;
    
    fs = data_struct.fsample;
    occipital_idx = find(ismember(data_struct.label, occipital_labels));
    
    feats = [];
    
    for t = 1:length(data_struct.trial)
        
        trial = data_struct.trial{t};
        signal = mean(trial(occipital_idx,:),1);
        
        [pxx, f] = pwelch(signal, [], [], [], fs);
        
        feat = zeros(1,nBands);
        
        for b = 1:nBands
            idx = f >= bands(b,1) & f <= bands(b,2);
            feat(b) = mean(pxx(idx));
        end
        
        feats = [feats; feat]; %#ok<AGROW>
    end
    
    feats = log(feats + 1e-8);
    feats = zscore(feats);
    
    features{c} = feats;
end

%% =========================
% COMPUTE CROSS-SIMILARITY
%% =========================

sim_matrix = zeros(3);

for i = 1:3
    for j = 1:3
        
        Xi = features{i};
        Xj = features{j};
        
        sims = [];
        
        for a = 1:size(Xi,1)
            for b = 1:size(Xj,1)
                sims = [sims; corr(Xi(a,:)', Xj(b,:)')]; %#ok<AGROW>
            end
        end
        
        sim_matrix(i,j) = mean(sims);
    end
end

%% =========================
% DISPLAY
%% =========================

disp('=== TRUE SIMILARITY MATRIX ===');
disp(sim_matrix);

figure;
imagesc(sim_matrix);
colorbar;
title('Trial-Level Similarity Matrix');
xticks(1:3); yticks(1:3);
xticklabels(labels_plot);
yticklabels(labels_plot);

%% =========================
% TRUE CII
%% =========================

CII = mean(sim_matrix(~eye(3)));

fprintf('\nTRUE CII: %.3f\n', CII);