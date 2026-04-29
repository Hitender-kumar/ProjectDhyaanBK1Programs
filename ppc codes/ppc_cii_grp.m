%% FINAL GROUP ANALYSIS (EXPERT vs CONTROL) - MULTI-BAND PPC
clear; clc; close all;

base_path = 'D:\ftavgdata';

%% =========================
% SUBJECT PAIRS + GENDER
%% =========================

pairs = {
    '019CKa','022SSP','M'
    '096MS','026HM','M'
    '040VS','100UK','M'
    '012GK','093AK','M'
    '095KM','075AD','M'
    '056PR','086AB','F'
    '052PR','082MS','F'
    '059MS','102AS','F'
    '013AR','064PK','F'
    '074KS','084AK','F'
};

nPairs = size(pairs,1);

%% =========================
% CHANNELS + BANDS
%% =========================

occipital_labels = {'O1','Oz','O2','POz','PO3','PO4','PO7','PO8','Iz'};
frontal_labels   = {'Fz','F1','F2','F3','F4'};

bands = [
    1 4;
    4 8;
    8 12;
    13 30
];

%% =========================
% MAIN LOOP
%% =========================

expert_base = []; expert_med = [];
control_base = []; control_med = [];

gender_labels = {};

for i = 1:nPairs

    exp_id = pairs{i,1};
    ctrl_id = pairs{i,2};
    gender = pairs{i,3};

    fprintf('\nPair %d: %s vs %s\n', i, exp_id, ctrl_id);

    [b1,m1] = compute_CII(fullfile(base_path,exp_id), occipital_labels, frontal_labels, bands);
    [b2,m2] = compute_CII(fullfile(base_path,ctrl_id), occipital_labels, frontal_labels, bands);

    expert_base = [expert_base; b1];
    expert_med  = [expert_med; m1];

    control_base = [control_base; b2];
    control_med  = [control_med; m2];

    gender_labels{i} = gender;
end

%% =========================
% GROUP RESULTS
%% =========================

fprintf('\n===== GROUP PPC MEANS =====\n');

fprintf('Expert Baseline: %.4f\n', mean(expert_base));
fprintf('Expert Meditation: %.4f\n', mean(expert_med));
fprintf('Control Baseline: %.4f\n', mean(control_base));
fprintf('Control Meditation: %.4f\n', mean(control_med));

%% =========================
% STATISTICS
%% =========================

[p_base,~] = ttest(expert_base, control_base);
[p_med,~]  = ttest(expert_med, control_med);

[p_exp_state,~] = ttest(expert_med, expert_base);
[p_ctrl_state,~] = ttest(control_med, control_base);

fprintf('\n===== P-VALUES =====\n');
fprintf('Baseline Expert vs Control: %.6e\n', p_base);
fprintf('Meditation Expert vs Control: %.6e\n', p_med);
fprintf('Expert Baseline vs Meditation: %.6e\n', p_exp_state);
fprintf('Control Baseline vs Meditation: %.6e\n', p_ctrl_state);

%% =========================
% EFFECT SIZE
%% =========================

d_base = mean(expert_base - control_base) / std(expert_base - control_base);
d_med  = mean(expert_med - control_med) / std(expert_med - control_med);

fprintf('\n===== EFFECT SIZE =====\n');
fprintf('Cohen d (Baseline): %.4f\n', d_base);
fprintf('Cohen d (Meditation): %.4f\n', d_med);

%% =========================
% INTERACTION
%% =========================

delta_exp = expert_med - expert_base;
delta_ctrl = control_med - control_base;

[p_interaction,~] = ttest(delta_exp, delta_ctrl);

fprintf('\n===== INTERACTION =====\n');
fprintf('Expert Improvement: %.4f\n', mean(delta_exp));
fprintf('Control Improvement: %.4f\n', mean(delta_ctrl));
fprintf('Interaction p-value: %.6e\n', p_interaction);

%% =========================
% GENDER ANALYSIS
%% =========================

male_idx = strcmp(gender_labels,'M');
female_idx = strcmp(gender_labels,'F');

fprintf('\n===== GENDER ANALYSIS =====\n');
fprintf('Male Expert Mean: %.4f\n', mean(expert_med(male_idx)));
fprintf('Female Expert Mean: %.4f\n', mean(expert_med(female_idx)));

%% =========================
% PLOTS
%% =========================

figure;
boxplot([expert_base, control_base, expert_med, control_med], ...
    'Labels', {'Exp-Base','Ctrl-Base','Exp-Med','Ctrl-Med'});
ylabel('PPC CII');
title('Expert vs Control (Multi-band PPC)');
grid on;

figure;
plot(expert_base - control_base, '-o'); hold on;
plot(expert_med - control_med, '-o');
legend('Baseline','Meditation');
title('Paired Differences (PPC)');
ylabel('CII Difference');
grid on;

figure;
boxplot([delta_exp, delta_ctrl], 'Labels', {'Expert','Control'});
ylabel('PPC Change');
title('Meditation Effect (PPC)');
grid on;

%% =========================
% FUNCTIONS
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

        % High-pass
        [b,a] = butter(3, 1/(fs/2), 'high');
        sig_occ = filtfilt(b,a,sig_occ);
        sig_fro = filtfilt(b,a,sig_fro);

        feat = zeros(1, size(bands,1));

        for bnd = 1:size(bands,1)

            band = bands(bnd,:);

            s1 = bandpass_filt(sig_fro, fs, band);
            s2 = bandpass_filt(sig_occ, fs, band);

            ph1 = angle(hilbert(s1));
            ph2 = angle(hilbert(s2));

            dphi = ph1 - ph2;

            feat(bnd) = compute_ppc(dphi);
        end

        feats = [feats; feat]; %#ok<AGROW>
    end

    feats = zscore(feats);
end

function sim = compute_similarity(X1,X2)

    vals = [];

    for i = 1:size(X1,1)
        for j = 1:size(X2,1)

            r = corr(X1(i,:)', X2(j,:)');
            vals = [vals; abs(r)]; %#ok<AGROW>

        end
    end

    sim = mean(vals);
end

function [CII_base, CII_med] = compute_CII(subj_path, occipital_labels, frontal_labels, bands)

    EO1 = extract_ppc_features(fullfile(subj_path,'EO1_ep_v8.mat'), occipital_labels, frontal_labels, bands);
    EC1 = extract_ppc_features(fullfile(subj_path,'EC1_ep_v8.mat'), occipital_labels, frontal_labels, bands);
    G1  = extract_ppc_features(fullfile(subj_path,'G1_ep_v8.mat'),  occipital_labels, frontal_labels, bands);

    EO2 = extract_ppc_features(fullfile(subj_path,'EO2_ep_v8.mat'), occipital_labels, frontal_labels, bands);
    EC2 = extract_ppc_features(fullfile(subj_path,'EC2_ep_v8.mat'), occipital_labels, frontal_labels, bands);
    G2  = extract_ppc_features(fullfile(subj_path,'G2_ep_v8.mat'),  occipital_labels, frontal_labels, bands);

    M1  = extract_ppc_features(fullfile(subj_path,'M1_ep_v8.mat'),  occipital_labels, frontal_labels, bands);
    M2  = extract_ppc_features(fullfile(subj_path,'M2_ep_v8.mat'),  occipital_labels, frontal_labels, bands);

    % Baseline
    b = [ ...
        compute_similarity(EO1,EC1), ...
        compute_similarity(EO1,G1), ...
        compute_similarity(EC1,G1)];

    CII_base = mean(b);

    % Meditation
    m = [ ...
        compute_similarity(EC2,EO2), ...
        compute_similarity(G2,EC2), ...
        compute_similarity(EO2,G2), ...
        compute_similarity(M2,EO2), ...
        compute_similarity(M2,EC2), ...
        compute_similarity(M2,G2)];

    CII_med = mean(m);
end

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