%% ================================================================
%  PAC_values_frontal_occipital_expert_vs_control.m

clear; clc;


data_root = 'D:\ftavgdata';  


pairs = {
    '019CKa', '022SSP', 'M'
    '096MS',  '026HM',  'M'
    '040VS',  '100UK',  'M'
    '012GK',  '093AK',  'M'
    '095KM',  '075AD',  'M'
    '056PR',  '086AB',  'F'
    '052PR',  '082MS',  'F'
    '059MS',  '102AS',  'F'
    '013AR',  '064PK',  'F'
    '074KS',  '084AK',  'F'
};

cond_files  = {'EO1_ep_v8.mat', 'EC1_ep_v8.mat', 'G1_ep_v8.mat', ...
               'M1_ep_v8.mat',  'G2_ep_v8.mat',  'EO2_ep_v8.mat', ...
               'EC2_ep_v8.mat', 'M2_ep_v8.mat'};
cond_labels = {'EO1','EC1','G1','M1','G2','EO2','EC2','M2'};

save_dir = fullfile(data_root, 'PAC_results_expert_vs_control');
if ~exist(save_dir, 'dir'), mkdir(save_dir); end

n_surr = 50;   
%% ===== FIXED PARAMETERS =====

fs     = 1000;
n_bins = 18;


alpha_lo = 8;  alpha_hi = 12;
gamma_lo = 30; gamma_hi = 80;
beta_lo  = 13; beta_hi  = 30;


occ_idx = [16,17,18, 46, 47, 48, 49, 50, 64];

% Fz=2, F3=3, F4=29, AFz=35, AF3=34, AF4=62, F1=36, F2=63
front_idx = [2,3,29,35,34,62,36,63];

nPairs      = size(pairs, 1);
nConditions = numel(cond_files);

expert_list  = pairs(:,1);
control_list = pairs(:,2);
sex_list     = pairs(:,3);


[b_hp, a_hp] = butter(3, 1/(fs/2),                    'high');
[b_ph, a_ph] = butter(3, [alpha_lo alpha_hi]/(fs/2),  'bandpass');
[b_ag, a_ag] = butter(3, [gamma_lo gamma_hi]/(fs/2),  'bandpass');
[b_ab, a_ab] = butter(3, [beta_lo  beta_hi ]/(fs/2),  'bandpass');

%% ===== GROUP STORAGE =====


exp_MI_1 = nan(nPairs, nConditions);
exp_MI_2 = nan(nPairs, nConditions);
exp_MI_3 = nan(nPairs, nConditions);
exp_MI_4 = nan(nPairs, nConditions);

ctrl_MI_1 = nan(nPairs, nConditions);
ctrl_MI_2 = nan(nPairs, nConditions);
ctrl_MI_3 = nan(nPairs, nConditions);
ctrl_MI_4 = nan(nPairs, nConditions);

exp_Z_1 = nan(nPairs, nConditions);
exp_Z_2 = nan(nPairs, nConditions);
exp_Z_3 = nan(nPairs, nConditions);
exp_Z_4 = nan(nPairs, nConditions);

ctrl_Z_1 = nan(nPairs, nConditions);
ctrl_Z_2 = nan(nPairs, nConditions);
ctrl_Z_3 = nan(nPairs, nConditions);
ctrl_Z_4 = nan(nPairs, nConditions);

%% ===== HELPER: process one subject =====

process_subject = @(subj_id) run_pac_subject(subj_id, data_root, ...
    cond_files, cond_labels, occ_idx, front_idx, ...
    b_hp, a_hp, b_ph, a_ph, b_ag, a_ag, b_ab, a_ab, ...
    n_bins, n_surr, save_dir);

%% ===== MAIN LOOP: pairs =====

for p = 1:nPairs

    exp_id  = expert_list{p};
    ctrl_id = control_list{p};
    sex     = sex_list{p};

    fprintf('\n===== Pair %d/%d : Expert=%s | Control=%s | Sex=%s =====\n', ...
            p, nPairs, exp_id, ctrl_id, sex);

    % ---- Process Expert ----
    fprintf('  >> Processing EXPERT: %s\n', exp_id);
    [m1,m2,m3,m4, z1,z2,z3,z4] = process_subject(exp_id);
    exp_MI_1(p,:) = m1;  exp_MI_2(p,:) = m2;
    exp_MI_3(p,:) = m3;  exp_MI_4(p,:) = m4;
    exp_Z_1(p,:)  = z1;  exp_Z_2(p,:)  = z2;
    exp_Z_3(p,:)  = z3;  exp_Z_4(p,:)  = z4;

    % ---- Process Control ----
    fprintf('  >> Processing CONTROL: %s\n', ctrl_id);
    [m1,m2,m3,m4, z1,z2,z3,z4] = process_subject(ctrl_id);
    ctrl_MI_1(p,:) = m1;  ctrl_MI_2(p,:) = m2;
    ctrl_MI_3(p,:) = m3;  ctrl_MI_4(p,:) = m4;
    ctrl_Z_1(p,:)  = z1;  ctrl_Z_2(p,:)  = z2;
    ctrl_Z_3(p,:)  = z3;  ctrl_Z_4(p,:)  = z4;

end % pairs

%% ===== SAVE GROUP SUMMARY =====

group_PAC.pairs       = pairs;
group_PAC.conditions  = cond_labels;

% Expert
group_PAC.exp_MI_frt_occ_gamma = exp_MI_1;
group_PAC.exp_MI_frt_occ_beta  = exp_MI_2;
group_PAC.exp_MI_occ_frt_gamma = exp_MI_3;
group_PAC.exp_MI_occ_frt_beta  = exp_MI_4;
group_PAC.exp_Z_frt_occ_gamma  = exp_Z_1;
group_PAC.exp_Z_frt_occ_beta   = exp_Z_2;
group_PAC.exp_Z_occ_frt_gamma  = exp_Z_3;
group_PAC.exp_Z_occ_frt_beta   = exp_Z_4;

% Control
group_PAC.ctrl_MI_frt_occ_gamma = ctrl_MI_1;
group_PAC.ctrl_MI_frt_occ_beta  = ctrl_MI_2;
group_PAC.ctrl_MI_occ_frt_gamma = ctrl_MI_3;
group_PAC.ctrl_MI_occ_frt_beta  = ctrl_MI_4;
group_PAC.ctrl_Z_frt_occ_gamma  = ctrl_Z_1;
group_PAC.ctrl_Z_frt_occ_beta   = ctrl_Z_2;
group_PAC.ctrl_Z_occ_frt_gamma  = ctrl_Z_3;
group_PAC.ctrl_Z_occ_frt_beta   = ctrl_Z_4;

save(fullfile(save_dir, 'PAC_frt_occ_expert_vs_control_summary.mat'), 'group_PAC');
fprintf('\nGroup summary saved.\n');

%% ===== STATISTICS: paired t-test + Wilcoxon + Cohen d =====

dir_names = {'Frt->Occ gamma (top-down)', 'Frt->Occ beta  (top-down)', ...
             'Occ->Frt gamma (bottom-up)','Occ->Frt beta  (bottom-up)'};

exp_all  = {exp_MI_1,  exp_MI_2,  exp_MI_3,  exp_MI_4};
ctrl_all = {ctrl_MI_1, ctrl_MI_2, ctrl_MI_3, ctrl_MI_4};

% Storage for stats results: [4 directions x nConditions]
stat_tval  = nan(4, nConditions);
stat_pval  = nan(4, nConditions);
stat_pWilc = nan(4, nConditions);
stat_cohd  = nan(4, nConditions);

fprintf('\n========== PAIRED STATISTICS (Expert vs Control) ==========\n');

for d = 1:4
    fprintf('\n--- %s ---\n', dir_names{d});
    fprintf('%-8s  %6s  %6s  %8s  %8s  %6s\n', ...
            'Cond','t-val','p(t)','p(Wilc)','Cohen-d','Sig');

    for c = 1:nConditions
        e_vec = exp_all{d}(:,c);
        k_vec = ctrl_all{d}(:,c);

        % Remove pairs where either is NaN
        valid_idx = ~isnan(e_vec) & ~isnan(k_vec);
        e_v = e_vec(valid_idx);
        k_v = k_vec(valid_idx);
        n_v = sum(valid_idx);

        if n_v < 3
            fprintf('%-8s  [insufficient data]\n', cond_labels{c});
            continue;
        end

        % Paired t-test
        [~, pt, ~, tstat_s] = ttest(e_v, k_v);
        tv = tstat_s.tstat;

        % Wilcoxon signed-rank
        [pw, ~] = signrank(e_v, k_v);

        % Cohen's d (paired: mean difference / SD of differences)
        diffs = e_v - k_v;
        cd    = mean(diffs) / (std(diffs) + eps);

        stat_tval(d,c)  = tv;
        stat_pval(d,c)  = pt;
        stat_pWilc(d,c) = pw;
        stat_cohd(d,c)  = cd;

        sig_str = '';
        if pt < 0.001, sig_str = '***';
        elseif pt < 0.01,  sig_str = '**';
        elseif pt < 0.05,  sig_str = '*';
        end

        fprintf('%-8s  %6.3f  %6.4f  %8.4f  %8.3f  %s\n', ...
                cond_labels{c}, tv, pt, pw, cd, sig_str);
    end
end

% Save stats
group_PAC.stat_tval  = stat_tval;
group_PAC.stat_pval  = stat_pval;
group_PAC.stat_pWilc = stat_pWilc;
group_PAC.stat_cohd  = stat_cohd;
save(fullfile(save_dir, 'PAC_frt_occ_expert_vs_control_summary.mat'), 'group_PAC');

%% ===== FIGURE 1: Group MI — Expert vs Control, 4 directions =====

col_exp  = [0.85 0.2  0.2 ];   % red   — Expert
col_ctrl = [0.2  0.4  0.75];   % blue  — Control

dir_short = {'Frt\alpha\rightarrowOcc\gamma (top-down)', ...
             'Frt\alpha\rightarrowOcc\beta  (top-down)', ...
             'Occ\alpha\rightarrowFrt\gamma (bottom-up)', ...
             'Occ\alpha\rightarrowFrt\beta  (bottom-up)'};

x_pos = 1:nConditions;
bar_w = 0.35;

figure('Name','Expert vs Control — Group MI','Position',[50 50 1300 560]);

for d = 1:4
    subplot(2,2,d);

    gm_e  = nanmean(exp_all{d},  1);
    sem_e = nanstd(exp_all{d},  0, 1) / sqrt(nPairs);
    gm_c  = nanmean(ctrl_all{d}, 1);
    sem_c = nanstd(ctrl_all{d}, 0, 1) / sqrt(nPairs);

    b1 = bar(x_pos - bar_w/2, gm_e,  bar_w, 'FaceColor', col_exp,  'EdgeColor','none'); hold on;
    b2 = bar(x_pos + bar_w/2, gm_c,  bar_w, 'FaceColor', col_ctrl, 'EdgeColor','none');

    errorbar(x_pos - bar_w/2, gm_e,  sem_e,  'k.','LineWidth',1.4,'CapSize',5);
    errorbar(x_pos + bar_w/2, gm_c,  sem_c,  'k.','LineWidth',1.4,'CapSize',5);

    % Mark significant conditions
    y_max = max([gm_e + sem_e, gm_c + sem_c]) * 1.15;
    for c = 1:nConditions
        if ~isnan(stat_pval(d,c)) && stat_pval(d,c) < 0.05
            sig_str = '*';
            if stat_pval(d,c) < 0.01,  sig_str = '**';  end
            if stat_pval(d,c) < 0.001, sig_str = '***'; end
            text(c, y_max, sig_str, 'HorizontalAlignment','center', ...
                 'FontSize',12,'Color','k','FontWeight','bold');
        end
    end

    xticks(x_pos); xticklabels(cond_labels);
    ylabel('KL Modulation Index');
    title(dir_short{d});
    legend([b1 b2], {'Expert','Control'}, 'Location','northwest','FontSize',7);
    grid on; box off;

    if d <= 2
        text(0.98,0.95,'TOP-DOWN','Units','normalized', ...
             'HorizontalAlignment','right','Color','r','FontWeight','bold','FontSize',8);
    else
        text(0.98,0.95,'BOTTOM-UP','Units','normalized', ...
             'HorizontalAlignment','right','Color','b','FontWeight','bold','FontSize',8);
    end
end

sgtitle('PAC Expert vs Control: Frontal \leftrightarrow Occipital  (* p<0.05, ** p<0.01, *** p<0.001)', ...
        'FontSize',12);

%% ===== FIGURE 2: Z-score heatmaps — Expert | Control side by side =====

if n_surr > 0

    exp_Z_all  = {exp_Z_1,  exp_Z_2,  exp_Z_3,  exp_Z_4};
    ctrl_Z_all = {ctrl_Z_1, ctrl_Z_2, ctrl_Z_3, ctrl_Z_4};

    all_Z = [cell2mat(exp_Z_all(:)); cell2mat(ctrl_Z_all(:))];
    z_max = max(all_Z(:), [], 'omitnan');

    figure('Name','Z-score Heatmap Expert vs Control','Position',[60 60 1400 500]);
    panel = 0;

    for d = 1:4
        for grp = 1:2   
            panel = panel + 1;
            subplot(4, 2, panel);

            if grp == 1
                dat     = exp_Z_all{d};
                grp_lbl = 'Expert';
                subj_ids = expert_list;
            else
                dat     = ctrl_Z_all{d};
                grp_lbl = 'Control';
                subj_ids = control_list;
            end

            imagesc(dat);
            colormap(hot);
            clim([0, max(z_max * 1.0, 2)]);
            colorbar;
            xticks(1:nConditions); xticklabels(cond_labels);
            yticks(1:nPairs);      yticklabels(subj_ids);
            xlabel('Condition'); ylabel('Subject');
            title(sprintf('%s | %s', grp_lbl, dir_short{d}), 'FontSize', 8);

            hold on;
            for r = 0.5:1:nPairs+0.5,      yline(r,'w-','LineWidth',0.5); end
            for cc = 0.5:1:nConditions+0.5, xline(cc,'w-','LineWidth',0.5); end
         
            for r = 1:nPairs
                for cc = 1:nConditions
                    if ~isnan(dat(r,cc)) && dat(r,cc) > 1.96
                        text(cc, r, '*','Color','c','FontSize',11, ...
                             'HorizontalAlignment','center', ...
                             'VerticalAlignment','middle','FontWeight','bold');
                    end
                end
            end
            hold off;
        end
    end

    sgtitle('PAC Z-scores: Expert (left) vs Control (right)  (* = Z>1.96)', ...
            'FontSize',11);
end

%% ===== FIGURE 3: Directionality Index — Expert vs Control =====

DI_gamma_exp  = (exp_MI_1  - exp_MI_3)  ./ (exp_MI_1  + exp_MI_3  + eps);
DI_beta_exp   = (exp_MI_2  - exp_MI_4)  ./ (exp_MI_2  + exp_MI_4  + eps);
DI_gamma_ctrl = (ctrl_MI_1 - ctrl_MI_3) ./ (ctrl_MI_1 + ctrl_MI_3 + eps);
DI_beta_ctrl  = (ctrl_MI_2 - ctrl_MI_4) ./ (ctrl_MI_2 + ctrl_MI_4 + eps);

figure('Name','Directionality Index Expert vs Control','Position',[100 100 1100 460]);

di_data   = {DI_gamma_exp, DI_gamma_ctrl; DI_beta_exp, DI_beta_ctrl};
di_titles = {'Gamma DI — Expert', 'Gamma DI — Control'; ...
             'Beta DI  — Expert',  'Beta DI  — Control'};
di_cols   = {[0.75 0.2 0.75], [0.3 0.6 0.9]; ...
             [0.75 0.2 0.75], [0.3 0.6 0.9]};

panel = 0;
for row = 1:2
    for grp = 1:2
        panel = panel + 1;
        subplot(2,2,panel);

        dat = di_data{row, grp};
        gm  = nanmean(dat, 1);
        sem = nanstd(dat, 0, 1) / sqrt(nPairs);

        bar(x_pos, gm, 0.6, 'FaceColor', di_cols{row,grp}, 'EdgeColor','none'); hold on;
        errorbar(x_pos, gm, sem, 'k.','LineWidth',1.5,'CapSize',6);
        for pp = 1:nPairs
            plot(x_pos, dat(pp,:), 'o-','Color',[0.75 0.75 0.75], ...
                 'MarkerSize',3,'LineWidth',0.8,'MarkerFaceColor',[0.75 0.75 0.75]);
        end
        yline(0,'k--','LineWidth',1.5);
        xticks(x_pos); xticklabels(cond_labels);
        ylabel('Directionality Index');
        title(di_titles{row,grp});
        ylim([-1 1]); grid on; box off;
    end
end

sgtitle({'Directionality Index: Expert vs Control', ...
         'Positive = Top-down dominant | Negative = Bottom-up dominant'}, ...
        'FontSize',11);

%% ===== FIGURE 4: Paired Difference Plots (Expert - Control) =====

figure('Name','Paired Differences Expert - Control','Position',[150 150 1200 500]);

for d = 1:4
    subplot(2,2,d);

    diff_mat = exp_all{d} - ctrl_all{d};   
    gm_diff  = nanmean(diff_mat, 1);
    sem_diff = nanstd(diff_mat, 0, 1) / sqrt(nPairs);

    
    bar_cols = repmat([0.6 0.6 0.6], nConditions, 1);
    bar_cols(gm_diff > 0, :) = repmat([0.85 0.2 0.2], sum(gm_diff > 0), 1);   
    bar_cols(gm_diff < 0, :) = repmat([0.2 0.4 0.75], sum(gm_diff < 0), 1);  

    for c = 1:nConditions
        bar(c, gm_diff(c), 0.6, 'FaceColor', bar_cols(c,:), 'EdgeColor','none'); hold on;
    end
    errorbar(x_pos, gm_diff, sem_diff, 'k.','LineWidth',1.4,'CapSize',5);

    
    for pp = 1:nPairs
        plot(x_pos, diff_mat(pp,:), 'o-', 'Color',[0.75 0.75 0.75], ...
             'MarkerSize',3,'LineWidth',0.8,'MarkerFaceColor',[0.75 0.75 0.75]);
    end

    yline(0,'k-','LineWidth',1.5);

    % Significance markers from paired t-test
    y_rng = max(abs([gm_diff + sem_diff, gm_diff - sem_diff]));
    for c = 1:nConditions
        if ~isnan(stat_pval(d,c)) && stat_pval(d,c) < 0.05
            sig_str = '*';
            if stat_pval(d,c) < 0.01,  sig_str = '**';  end
            if stat_pval(d,c) < 0.001, sig_str = '***'; end
            text(c, y_rng * 1.15, sig_str, 'HorizontalAlignment','center', ...
                 'FontSize',12,'Color','k','FontWeight','bold');
        end
    end

    xticks(x_pos); xticklabels(cond_labels);
    ylabel('\Delta MI (Expert - Control)');
    title(dir_short{d});
    grid on; box off;
    text(0.02, 0.93,'Expert > Control','Units','normalized','Color','r','FontSize',8);
    text(0.02, 0.07,'Control > Expert','Units','normalized','Color','b','FontSize',8);
end

sgtitle({'Paired Difference: Expert - Control MI', ...
         '(* p<0.05, ** p<0.01, *** p<0.001, paired t-test)'}, 'FontSize',11);

fprintf('\nAll done. Results saved to: %s\n', save_dir);

%% ================================================================


function [mean1, mean2, mean3, mean4, Z1, Z2, Z3, Z4] = ...
    run_pac_subject(subj_id, data_root, cond_files, cond_labels, ...
                    occ_idx, front_idx, ...
                    b_hp, a_hp, b_ph, a_ph, b_ag, a_ag, b_ab, a_ab, ...
                    n_bins, n_surr, save_dir)

    nConditions = numel(cond_files);
    alpha_lo = 8; alpha_hi = 12;
    gamma_lo = 30; gamma_hi = 80;
    beta_lo  = 13; beta_hi  = 30;
    fs = 1000;

    subj_dir = fullfile(data_root, subj_id);

    mean1 = nan(1,nConditions);  mean2 = nan(1,nConditions);
    mean3 = nan(1,nConditions);  mean4 = nan(1,nConditions);
    Z1    = nan(1,nConditions);  Z2    = nan(1,nConditions);
    Z3    = nan(1,nConditions);  Z4    = nan(1,nConditions);

    MI_cell_1 = cell(nConditions,1);  MI_cell_2 = cell(nConditions,1);
    MI_cell_3 = cell(nConditions,1);  MI_cell_4 = cell(nConditions,1);

    for c = 1:nConditions

        fpath = fullfile(subj_dir, cond_files{c});
        if ~exist(fpath, 'file')
            fprintf('    [%s | %s] File not found — skipping\n', subj_id, cond_labels{c});
            continue;
        end

        raw = load(fpath);
        ft  = raw.data;

        trials_cell = ft.trial;
        nTrials     = length(trials_cell);
        bad_elecs   = ft.badElecs(:);

        occ_clean   = setdiff(occ_idx,   bad_elecs);
        front_clean = setdiff(front_idx, bad_elecs);

        if isempty(occ_clean)
            fprintf('    [%s | %s] All occ channels bad — skipping\n', subj_id, cond_labels{c});
            continue;
        end
        if isempty(front_clean)
            fprintf('    [%s | %s] All frontal channels bad — skipping\n', subj_id, cond_labels{c});
            continue;
        end

        fprintf('    [%s | %s] %d trials | occ=%d ch | frt=%d ch\n', ...
                subj_id, cond_labels{c}, nTrials, numel(occ_clean), numel(front_clean));

        MI_1 = zeros(1,nTrials);  MI_2 = zeros(1,nTrials);
        MI_3 = zeros(1,nTrials);  MI_4 = zeros(1,nTrials);

        n_total_surr = nTrials * n_surr;
        surr_1 = zeros(1,n_total_surr);  surr_2 = zeros(1,n_total_surr);
        surr_3 = zeros(1,n_total_surr);  surr_4 = zeros(1,n_total_surr);
        sc = 0;

        for t = 1:nTrials

            trial = ft.trial{t};
            n_t   = size(trial, 2);

            occ_sig = filtfilt(b_hp, a_hp, mean(trial(occ_clean,   :), 1));
            frt_sig = filtfilt(b_hp, a_hp, mean(trial(front_clean, :), 1));

            ph_frt = angle(hilbert(filtfilt(b_ph, a_ph, frt_sig)));
            ph_occ = angle(hilbert(filtfilt(b_ph, a_ph, occ_sig)));

            env_occ_g = abs(hilbert(filtfilt(b_ag, a_ag, occ_sig)));
            env_occ_b = abs(hilbert(filtfilt(b_ab, a_ab, occ_sig)));
            env_frt_g = abs(hilbert(filtfilt(b_ag, a_ag, frt_sig)));
            env_frt_b = abs(hilbert(filtfilt(b_ab, a_ab, frt_sig)));

            MI_1(t) = tort_MI(ph_frt, env_occ_g, n_bins);
            MI_2(t) = tort_MI(ph_frt, env_occ_b, n_bins);
            MI_3(t) = tort_MI(ph_occ, env_frt_g, n_bins);
            MI_4(t) = tort_MI(ph_occ, env_frt_b, n_bins);

            if n_surr > 0
                for sv = 1:n_surr
                    shift = randi([round(n_t * 0.1), round(n_t * 0.9)]);
                    sc = sc + 1;
                    surr_1(sc) = tort_MI(ph_frt, circshift(env_occ_g, shift), n_bins);
                    surr_2(sc) = tort_MI(ph_frt, circshift(env_occ_b, shift), n_bins);
                    surr_3(sc) = tort_MI(ph_occ, circshift(env_frt_g, shift), n_bins);
                    surr_4(sc) = tort_MI(ph_occ, circshift(env_frt_b, shift), n_bins);
                end
            end

        end % trials

        mean1(c) = mean(MI_1);  mean2(c) = mean(MI_2);
        mean3(c) = mean(MI_3);  mean4(c) = mean(MI_4);

        MI_cell_1{c} = MI_1;  MI_cell_2{c} = MI_2;
        MI_cell_3{c} = MI_3;  MI_cell_4{c} = MI_4;

        if n_surr > 0 && sc > 0
            s1 = surr_1(1:sc);  s2 = surr_2(1:sc);
            s3 = surr_3(1:sc);  s4 = surr_4(1:sc);

            Z1(c) = (mean1(c) - mean(s1)) / (std(s1) + eps);
            Z2(c) = (mean2(c) - mean(s2)) / (std(s2) + eps);
            Z3(c) = (mean3(c) - mean(s3)) / (std(s3) + eps);
            Z4(c) = (mean4(c) - mean(s4)) / (std(s4) + eps);
        end

        fprintf('      Frt->Occ gamma: MI=%.5f  Z=%.2f\n', mean1(c), Z1(c));
        fprintf('      Frt->Occ beta : MI=%.5f  Z=%.2f\n', mean2(c), Z2(c));
        fprintf('      Occ->Frt gamma: MI=%.5f  Z=%.2f\n', mean3(c), Z3(c));
        fprintf('      Occ->Frt beta : MI=%.5f  Z=%.2f\n', mean4(c), Z4(c));

    end % conditions

    % Per-subject save
    pac_subj.subject_id        = subj_id;
    pac_subj.conditions        = cond_labels;
    pac_subj.MI_frt_occ_gamma  = MI_cell_1;
    pac_subj.MI_frt_occ_beta   = MI_cell_2;
    pac_subj.MI_occ_frt_gamma  = MI_cell_3;
    pac_subj.MI_occ_frt_beta   = MI_cell_4;
    pac_subj.mean_frt_occ_gamma = mean1;
    pac_subj.mean_frt_occ_beta  = mean2;
    pac_subj.mean_occ_frt_gamma = mean3;
    pac_subj.mean_occ_frt_beta  = mean4;
    pac_subj.Z_frt_occ_gamma    = Z1;
    pac_subj.Z_frt_occ_beta     = Z2;
    pac_subj.Z_occ_frt_gamma    = Z3;
    pac_subj.Z_occ_frt_beta     = Z4;
    pac_subj.occ_idx    = occ_idx;
    pac_subj.front_idx  = front_idx;
    pac_subj.freq_bands = struct('alpha',[8 12], 'gamma',[30 80], 'beta',[13 30]);

    save(fullfile(save_dir, sprintf('PAC_frt_occ_%s.mat', subj_id)), 'pac_subj');
    fprintf('    Saved PAC_frt_occ_%s.mat\n', subj_id);

end 
%% ===== LOCAL FUNCTION: tort_MI =====

function MI = tort_MI(phase, amp, n_bins)
    edges      = linspace(-pi, pi, n_bins + 1);
    amp_binned = zeros(1, n_bins);
    for k = 1:n_bins
        idx = phase >= edges(k) & phase < edges(k+1);
        if any(idx)
            amp_binned(k) = mean(amp(idx));
        end
    end
    s        = sum(amp_binned) + eps;
    amp_dist = amp_binned / s;
    H        = -sum(amp_dist .* log(amp_dist + eps));
    MI       = (log(n_bins) - H) / log(n_bins);
end