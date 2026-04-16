%% PHASE 3: FULL 8-CONDITION + MEDITATION ANALYSIS (13AR)
clear; clc; close all;

%% =========================
% PATH (EDIT IF NEEDED)
%% =========================

base_path = 'C:\Users\VAIBHAV\Desktop\NSP_GrantProposal\ProjectDhyaanBK1Programs-master\data\ftDataAvgRef\013AR\';

files = {
    'EO1_ep_v8.mat'   % 1
    'EC1_ep_v8.mat'   % 2
    'G1_ep_v8.mat'    % 3
    'M1_ep_v8.mat'    % 4
    'G2_ep_v8.mat'    % 5
    'EO2_ep_v8.mat'   % 6
    'EC2_ep_v8.mat'   % 7
    'M2_ep_v8.mat'    % 8
};

labels = {'EO1','EC1','G1','M1','G2','EO2','EC2','M2'};

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

%% =========================
% FEATURE EXTRACTION FUNCTION
%% =========================

function feats = extract_features(file_name, occipital_labels, bands)

    data = load(file_name);
    data_struct = data.data;

    fs = data_struct.fsample;
    occipital_idx = find(ismember(data_struct.label, occipital_labels));

    feats = [];

    for t = 1:length(data_struct.trial)

        trial = data_struct.trial{t};
        signal = mean(trial(occipital_idx,:),1);

        [pxx, f] = pwelch(signal, [], [], [], fs);

        feat = zeros(1, size(bands,1));

        for b = 1:size(bands,1)
            idx = f >= bands(b,1) & f <= bands(b,2);
            feat(b) = mean(pxx(idx));
        end

        feats = [feats; feat]; %#ok<AGROW>
    end

    feats = log(feats + 1e-8);
    feats = zscore(feats);
end

%% =========================
% EXTRACT FEATURES
%% =========================

features = cell(length(files),1);

for i = 1:length(files)
    fprintf('Processing %s\n', files{i});
    features{i} = extract_features(fullfile(base_path, files{i}), occipital_labels, bands);
end

%% =========================
% SIMILARITY FUNCTION
%% =========================

function sim = compute_similarity(X1, X2)

    sims = [];

    for i = 1:size(X1,1)
        for j = 1:size(X2,1)
            r = corr(X1(i,:)', X2(j,:)');
            sims = [sims; abs(r)]; %#ok<AGROW>
        end
    end

    sim = mean(sims);
end

%% =========================
% FULL SIMILARITY MATRIX
%% =========================

n = length(files);
sim_matrix = zeros(n);

for i = 1:n
    for j = 1:n
        sim_matrix(i,j) = compute_similarity(features{i}, features{j});
    end
end

%% =========================
% VISUALIZE
%% =========================

figure;
imagesc(sim_matrix);
colorbar;
title('8-Condition Similarity Matrix (Full-band)');
xticks(1:n); yticks(1:n);
xticklabels(labels);
yticklabels(labels);

%% =========================
% MEDITATION ANALYSIS
%% =========================

% Indices
EO1 = 1; EC1 = 2; G1 = 3; M1 = 4;
G2  = 5; EO2 = 6; EC2 = 7; M2 = 8;

%% --- Compare M1 ---
sim_M1_EO = sim_matrix(M1, EO1);
sim_M1_G  = sim_matrix(M1, G1);
sim_M1_EC = sim_matrix(M1, EC1);

%% --- Compare M2 ---
sim_M2_EO = sim_matrix(M2, EO2);
sim_M2_G  = sim_matrix(M2, G2);
sim_M2_EC = sim_matrix(M2, EC2);

fprintf('\n===== MEDITATION SIMILARITY =====\n');

fprintf('\n--- M1 ---\n');
fprintf('M1 vs EO1: %.4f\n', sim_M1_EO);
fprintf('M1 vs G1 : %.4f\n', sim_M1_G);
fprintf('M1 vs EC1: %.4f\n', sim_M1_EC);

fprintf('\n--- M2 ---\n');
fprintf('M2 vs EO2: %.4f\n', sim_M2_EO);
fprintf('M2 vs G2 : %.4f\n', sim_M2_G);
fprintf('M2 vs EC2: %.4f\n', sim_M2_EC);

%% =========================
% GLOBAL CII COMPARISON
%% =========================

baseline_pairs = [
    sim_matrix(1,2), % EO1-EC1
    sim_matrix(1,3), % EO1-G1
    sim_matrix(2,3)  % EC1-G1
];

meditation_pairs = [
    sim_M1_EO, sim_M1_G, sim_M1_EC, ...
    sim_M2_EO, sim_M2_G, sim_M2_EC
];

CII_baseline = mean(baseline_pairs);
CII_meditation = mean(meditation_pairs);

fprintf('\n===== CII COMPARISON =====\n');
fprintf('Baseline CII: %.4f\n', CII_baseline);
fprintf('Meditation CII: %.4f\n', CII_meditation);

%% =========================
% BAR PLOT
%% =========================

figure;
bar([CII_baseline, CII_meditation]);
set(gca,'XTickLabel',{'Baseline','Meditation'});
ylabel('Similarity');
title('CII: Baseline vs Meditation');
grid on;