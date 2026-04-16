%% PHASE 2: TEMPORAL + ALPHA + FULL-BAND CII (FIXED) - 13AR
clear; clc; close all;

%% =========================
% PATH (EDIT IF NEEDED)
%% =========================

base_path = 'C:\Users\VAIBHAV\Desktop\NSP_GrantProposal\ProjectDhyaanBK1Programs-master\data\ftDataAvgRef\013AR\';

files = {
    'EC1_ep_v8.mat'
    'EC2_ep_v8.mat'
    'EO1_ep_v8.mat'
    'EO2_ep_v8.mat'
    'G1_ep_v8.mat'
    'G2_ep_v8.mat'
};

labels = {'EC1','EC2','EO1','EO2','G1','G2'};

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

alpha_band = [8 12];

%% =========================
% FEATURE EXTRACTION FUNCTION
%% =========================

function [full_feats, alpha_feats] = extract_features(file_name, occipital_labels, bands, alpha_band)

    data = load(file_name);
    data_struct = data.data;

    fs = data_struct.fsample;
    occipital_idx = find(ismember(data_struct.label, occipital_labels));

    full_feats = [];
    alpha_feats = [];

    for t = 1:length(data_struct.trial)

        trial = data_struct.trial{t};
        signal = mean(trial(occipital_idx,:),1);

        [pxx, f] = pwelch(signal, [], [], [], fs);

        % FULL BAND FEATURES
        feat = zeros(1, size(bands,1));
        for b = 1:size(bands,1)
            idx = f >= bands(b,1) & f <= bands(b,2);
            feat(b) = mean(pxx(idx));
        end

        % ALPHA ONLY
        idx_alpha = f >= alpha_band(1) & f <= alpha_band(2);
        alpha_val = mean(pxx(idx_alpha));

        full_feats = [full_feats; feat]; %#ok<AGROW>
        alpha_feats = [alpha_feats; alpha_val]; %#ok<AGROW>
    end

    % Normalize
    full_feats = log(full_feats + 1e-8);
    full_feats = zscore(full_feats);

    alpha_feats = log(alpha_feats + 1e-8);
    alpha_feats = zscore(alpha_feats);
end

%% =========================
% EXTRACT FEATURES
%% =========================

full_features = cell(6,1);
alpha_features = cell(6,1);

for i = 1:6
    fprintf('Processing %s\n', files{i});
    [full_features{i}, alpha_features{i}] = extract_features( ...
        fullfile(base_path, files{i}), occipital_labels, bands, alpha_band);
end

%% =========================
% SIMILARITY FUNCTION (FIXED)
%% =========================

function sim = compute_similarity(X1, X2)

    sims = [];

    for i = 1:size(X1,1)
        for j = 1:size(X2,1)
            
            if size(X1,2) == 1
                % ALPHA CASE (1D)
                d = abs(X1(i) - X2(j));
                s = exp(-d);
            else
                % FULL-BAND CASE (multi-D)
                r = corr(X1(i,:)', X2(j,:)');
                s = abs(r);
            end
            
            sims = [sims; s]; %#ok<AGROW>
        end
    end

    sim = mean(sims);
end

%% =========================
% TEMPORAL STABILITY (FULL-BAND)
%% =========================

sim_EC = compute_similarity(full_features{1}, full_features{2});
sim_EO = compute_similarity(full_features{3}, full_features{4});
sim_G  = compute_similarity(full_features{5}, full_features{6});

fprintf('\n===== TEMPORAL STABILITY (FULL-BAND) =====\n');
fprintf('EC1 vs EC2: %.4f\n', sim_EC);
fprintf('EO1 vs EO2: %.4f\n', sim_EO);
fprintf('G1 vs G2:  %.4f\n', sim_G);

%% =========================
% ALPHA CII
%% =========================

alpha_EC = compute_similarity(alpha_features{1}, alpha_features{2});
alpha_EO = compute_similarity(alpha_features{3}, alpha_features{4});
alpha_G  = compute_similarity(alpha_features{5}, alpha_features{6});

fprintf('\n===== ALPHA CII =====\n');
fprintf('EC (Alpha): %.4f\n', alpha_EC);
fprintf('EO (Alpha): %.4f\n', alpha_EO);
fprintf('G  (Alpha): %.4f\n', alpha_G);

%% =========================
% FULL-BAND CII
%% =========================

full_EC = sim_EC;
full_EO = sim_EO;
full_G  = sim_G;

fprintf('\n===== FULL-BAND CII =====\n');
fprintf('EC (Full): %.4f\n', full_EC);
fprintf('EO (Full): %.4f\n', full_EO);
fprintf('G  (Full): %.4f\n', full_G);

%% =========================
% COMPARISON PLOT
%% =========================

figure;

subplot(1,2,1)
bar([alpha_EC, alpha_EO, alpha_G]);
title('Alpha CII');
set(gca,'XTickLabel',{'EC','EO','G'});
ylabel('Similarity');
grid on;

subplot(1,2,2)
bar([full_EC, full_EO, full_G]);
title('Full-band CII');
set(gca,'XTickLabel',{'EC','EO','G'});
ylabel('Similarity');
grid on;

%% =========================
% FINAL SUMMARY
%% =========================

fprintf('\n===== SUMMARY =====\n');

fprintf('Alpha vs Full-band comparison:\n');
fprintf('EC: Alpha=%.4f | Full=%.4f\n', alpha_EC, full_EC);
fprintf('EO: Alpha=%.4f | Full=%.4f\n', alpha_EO, full_EO);
fprintf('G : Alpha=%.4f | Full=%.4f\n', alpha_G, full_G);