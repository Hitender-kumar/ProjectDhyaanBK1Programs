%% FULL 8-CONDITION + PPC CONNECTIVITY ANALYSIS (MULTI-BAND)
clear; clc; close all;

%% =========================
% PATH
%% =========================

base_path = 'D:\ftavgdata\013AR';

files = {
    'EO1_ep_v8.mat'
    'EC1_ep_v8.mat'
    'G1_ep_v8.mat'
    'M1_ep_v8.mat'
    'G2_ep_v8.mat'
    'EO2_ep_v8.mat'
    'EC2_ep_v8.mat'
    'M2_ep_v8.mat'
};

labels = {'EO1','EC1','G1','M1','G2','EO2','EC2','M2'};

%% =========================
% CHANNELS
%% =========================

occipital_labels = {'O1','Oz','O2','POz','PO3','PO4','PO7','PO8','Iz'};
frontal_labels   = {'Fz','F1','F2','F3','F4'};

%% =========================
% FREQUENCY BANDS
%% =========================

bands = [
    1 4;
    4 8;
    8 12;
    13 30
];

%% =========================
% FEATURE EXTRACTION (MULTI-BAND PPC)
%% =========================

function feats = extract_ppc_features(file_name, occipital_labels, frontal_labels, bands)

    data = load(file_name);
    data_struct = data.data;

    fs = data_struct.fsample;

    occ_idx = find(ismember(data_struct.label, occipital_labels));
    fro_idx = find(ismember(data_struct.label, frontal_labels));

    feats = [];

    for t = 1:length(data_struct.trial)

        trial = data_struct.trial{t};

        sig_occ = mean(trial(occ_idx,:),1);
        sig_fro = mean(trial(fro_idx,:),1);

        % High-pass filter
        [b,a] = butter(3, 1/(fs/2), 'high');
        sig_occ = filtfilt(b,a,sig_occ);
        sig_fro = filtfilt(b,a,sig_fro);

        feat = zeros(1, size(bands,1));

        for bnd = 1:size(bands,1)

            band = bands(bnd,:);

            % Bandpass filter
            s1 = bandpass_filt(sig_fro, fs, band);
            s2 = bandpass_filt(sig_occ, fs, band);

            % Phase extraction
            ph1 = angle(hilbert(s1));
            ph2 = angle(hilbert(s2));

            dphi = ph1 - ph2;

            % PPC computation
            feat(bnd) = compute_ppc(dphi);
        end

        feats = [feats; feat]
    end

    % Normalize
    feats = zscore(feats);
end

%% =========================
% EXTRACT FEATURES
%% =========================

features = cell(length(files),1);

for i = 1:length(files)
    fprintf('Processing %s\n', files{i});
    features{i} = extract_ppc_features( ...
        fullfile(base_path, files{i}), ...
        occipital_labels, frontal_labels, bands);
end

%% =========================
% SIMILARITY FUNCTION (MULTI-D)
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
title('8-Condition Similarity Matrix (Multi-band PPC)');
xticks(1:n); yticks(1:n);
xticklabels(labels);
yticklabels(labels);

%% =========================
% MEDITATION ANALYSIS
%% =========================

EO1 = 1; EC1 = 2; G1 = 3; M1 = 4;
G2  = 5; EO2 = 6; EC2 = 7; M2 = 8;

sim_M1_EO = sim_matrix(M1, EO1);
sim_M1_G  = sim_matrix(M1, G1);
sim_M1_EC = sim_matrix(M1, EC1);

sim_M2_EO = sim_matrix(M2, EO2);
sim_M2_G  = sim_matrix(M2, G2);
sim_M2_EC = sim_matrix(M2, EC2);

fprintf('\n===== PPC MEDITATION SIMILARITY =====\n');

fprintf('\n--- M1 ---\n');
fprintf('M1 vs EO1: %.4f\n', sim_M1_EO);
fprintf('M1 vs G1 : %.4f\n', sim_M1_G);
fprintf('M1 vs EC1: %.4f\n', sim_M1_EC);

fprintf('\n--- M2 ---\n');
fprintf('M2 vs EO2: %.4f\n', sim_M2_EO);
fprintf('M2 vs G2 : %.4f\n', sim_M2_G);
fprintf('M2 vs EC2: %.4f\n', sim_M2_EC);

%% =========================
% GLOBAL CII
%% =========================

baseline_pairs = [
    sim_matrix(1,2)
    sim_matrix(1,3)
    sim_matrix(2,3)
];

meditation_pairs = [
    sim_M1_EO, sim_M1_G, sim_M1_EC, ...
    sim_M2_EO, sim_M2_G, sim_M2_EC
];

CII_baseline = mean(baseline_pairs);
CII_meditation = mean(meditation_pairs);

fprintf('\n===== PPC CII =====\n');
fprintf('Baseline: %.4f\n', CII_baseline);
fprintf('Meditation: %.4f\n', CII_meditation);

%% =========================
% BAR PLOT
%% =========================

figure;
bar([CII_baseline, CII_meditation]);
set(gca,'XTickLabel',{'Baseline','Meditation'});
ylabel('Similarity');
title('PPC CII (Multi-band)');
grid on;

%% =========================
% FILTER + PPC FUNCTIONS
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