% =========================================================================
% NSWOA 多目标姿态角优化
% 用于求解31个刀位点最优 (alpha_i, beta_i)，同时最小化 F1 和 F2
% 参考文献格式输出：Pareto前沿、收敛过程、最优折衷解细节
% =========================================================================
clear; clc; close all;

%% ---------- 常数与问题参数设置 ----------
RT = 4;                         % 常数 R_T
c0 = 1;                         % 常数 c_0
N_points = 31;                  % 刀位点数
dim = 2 * N_points;             % 决策变量维数 (31个alpha + 31个beta)

% 点数据: [xi, yi, zi, theta_i, alpha_min, alpha_max, beta_min, beta_max]
points_data = [
    20, -21.99999,   0,        0,          41.40962211, 138.5903779,   8.565282944, 48.59037789;
    20, -21.8795,   -2.2996,   5.999927178, 41.40962211, 138.5903779,  11.38537583,  54.59030507;
    20, -21.5192,   -4.5741,  12.0001346,  41.40962211, 138.5903779,  14.22831014,  60.59051249;
    20, -20.9232,   -6.7984,  18.0000996,  41.40962211, 138.5903779,  17.09324644,  66.59047749;
    20, -20.098,    -8.9482,  23.99998545, 41.40962211, 138.5903779,  19.97999322,  72.59036334;
    20, -19.0526,  -11,       29.99994646, 41.40962211, 138.5357316,  22.88865588,  78.59032435;
    20, -17.7984,  -12.9313,  36.00001152, 41.40962211, 134.9958803,  25.81958871,  84.59038941;
    20, -16.3492,  -14.7209,  42.00002748, 41.40962211, 132.0298962,  28.77341238,  90.59040537;
    20, -14.7209,  -16.3492,  47.99997252, 41.40962211, 129.5748898,  31.75120048,  96.59035041;
    20, -12.9313,  -17.7984,  53.99998848, 41.40962211, 127.5699073,  34.75455286, 102.5903664;
    20, -11,       -19.0526,  60.00005354, 41.40962211, 125.9609898,  37.78551104, 108.5904314;
    20,  -8.9482,  -20.098,   66.00001455, 41.40962211, 124.7030814,  40.84667876, 114.5903924;
    20,  -6.7984,  -20.9232,  71.9999004,  41.40962211, 123.7602199,  43.94159802, 120.5902783;
    20,  -4.5741,  -21.5192,  77.9998654,  41.40962211, 123.1053379,  47.0748956,  123.8050332;
    20,  -2.2996,  -21.8795,  84.00007282, 41.40962211, 122.7196641,  46.9273849,  126.9299426;
    20,   0,       -22,       90,          41.40962211, 122.5923577,  49.98297374, 130.0170263;
    20,   2.2996,  -21.8795,  95.99992718, 41.40962211, 122.7196641,  53.07005737, 133.0726151;
    20,   4.5741,  -21.5192, 102.0001346,  41.40962211, 123.1053379,  56.19496681, 132.9251044;
    20,   6.7984,  -20.9232, 108.0000996,  41.40962211, 123.7602199,  59.40972171, 136.058402;
    20,   8.9482,  -20.098,  113.9999854,  41.40962211, 124.7030814,  65.40960755, 139.1533212;
    20,  11,       -19.0526, 119.9999465,  41.40962211, 125.9609898,  71.40956857, 142.214489;
    20,  12.9313,  -17.7984, 126.0000115,  41.40962211, 127.5699073,  77.40963363, 145.2454471;
    20,  14.7209,  -16.3492, 132.0000275,  41.40962211, 129.5748898,  83.40964959, 148.2487995;
    20,  16.3492,  -14.7209, 137.9999725,  41.40962211, 132.0298962,  89.40959463, 151.2265876;
    20,  17.7984,  -12.9313, 143.9999885,  41.40962211, 134.9958803,  95.40961059, 154.1804113;
    20,  19.0526,  -11,      150.0000535,  41.40962211, 138.5357316, 101.4096757,  157.1113441;
    20,  20.098,   -8.9482,  156.0000146,  41.40962211, 138.5903779, 107.4096367,  160.0200068;
    20,  20.9232,  -6.7984,  161.9999004,  41.40962211, 138.5903779, 113.4095225,  162.9067536;
    20,  21.5192,  -4.5741,  167.9998654,  41.40962211, 138.5903779, 119.4094875,  165.7716899;
    20,  21.8795,  -2.2996,  174.0000728,  41.40962211, 138.5903779, 125.4096949,  168.6146242;
    20,  21.99999,   0,      180,         41.40962211, 138.5903779, 131.4096221,  171.4347171
];

% 提取 theta (度)
theta_deg = points_data(:,4)';

% 决策变量上下界 (度)
lb_alpha = points_data(:,5)';   % alpha 下界
ub_alpha = points_data(:,6)';   % alpha 上界
lb_beta  = points_data(:,7)';   % beta  下界
ub_beta  = points_data(:,8)';   % beta  上界

lb = [lb_alpha, lb_beta];       % 1x62 下界
ub = [ub_alpha, ub_beta];       % 1x62 上界

%% ---------- NSWOA 算法参数设置 ----------
pop_size = 100;                 % 种群规模
max_iter = 200;                 % 最大迭代次数
n_obj = 2;                      % 目标函数个数

% WOA 参数
b_woa = 1;                      % 螺旋形状常数

% 用于收敛过程记录
pareto_front_history = cell(max_iter, 1);  % 每代第一前沿的目标值

%% ---------- 运行 NSWOA 算法 ----------
fprintf('NSWOA 优化开始...\n');
[final_pop, final_obj, pareto_front_history] = NSWOA(@obj_fun, dim, lb, ub, ...
    pop_size, max_iter, n_obj, b_woa, theta_deg, RT, c0);

%% ---------- 结果提取与可视化 ----------
% 最终种群的非支配排序，提取 Pareto 前沿
[front_rank, crowd_dist] = NonDominatedSorting(final_obj);
pareto_idx = (front_rank == 1);
pareto_solutions = final_pop(pareto_idx, :);
pareto_objectives = final_obj(pareto_idx, :);

% 输出前沿信息
fprintf('\n========== 最终 Pareto 前沿 ==========\n');
fprintf('前沿解个数: %d\n', size(pareto_objectives,1));
fprintf('F1 范围: [%.6f, %.6f]\n', min(pareto_objectives(:,1)), max(pareto_objectives(:,1)));
fprintf('F2 范围: [%.6f, %.6f]\n', min(pareto_objectives(:,2)), max(pareto_objectives(:,2)));

%% 图1: 最终 Pareto 前沿
figure('Color','white','Position',[100,100,800,600]);
plot(pareto_objectives(:,1), pareto_objectives(:,2), 'ro', 'MarkerSize', 8, ...
    'MarkerFaceColor','r');
xlabel('F_1 (负平均 I)','FontSize',12);
ylabel('F_2 (平均 \gamma^2)','FontSize',12);
title('NSWOA 获得的 Pareto 前沿','FontSize',14);
grid on; set(gca,'FontSize',11);
saveas(gcf, 'Pareto_Front.png');

%% 图2: 进化过程 Pareto 前沿变化 (选取第1代、中代、最后一代)
figure('Color','white','Position',[100,100,800,600]);
generations_to_plot = [1, round(max_iter/2), max_iter];
colors = {'b','g','r'};
legend_str = {};
hold on;
for k = 1:length(generations_to_plot)
    gen = generations_to_plot(k);
    pf = pareto_front_history{gen};
    if ~isempty(pf)
        plot(pf(:,1), pf(:,2), 'o', 'Color', colors{k}, ...
            'MarkerSize', 6, 'MarkerFaceColor', colors{k});
        legend_str{end+1} = sprintf('第 %d 代', gen);
    end
end
hold off;
xlabel('F_1','FontSize',12); ylabel('F_2','FontSize',12);
title('不同迭代阶段的 Pareto 前沿','FontSize',14);
legend(legend_str,'Location','best');
grid on; set(gca,'FontSize',11);
saveas(gcf, 'Pareto_Evolution.png');

%% 图3: 收敛曲线 (每代第一前沿的 F1 和 F2 平均值)
avg_F1 = zeros(max_iter,1);
avg_F2 = zeros(max_iter,1);
for gen = 1:max_iter
    pf = pareto_front_history{gen};
    if ~isempty(pf)
        avg_F1(gen) = mean(pf(:,1));
        avg_F2(gen) = mean(pf(:,2));
    else
        avg_F1(gen) = NaN;
        avg_F2(gen) = NaN;
    end
end

figure('Color','white','Position',[100,100,800,400]);
subplot(1,2,1);
plot(1:max_iter, avg_F1, 'b-','LineWidth',1.5);
xlabel('迭代次数'); ylabel('第一前沿平均 F_1');
title('F_1 收敛曲线'); grid on;
subplot(1,2,2);
plot(1:max_iter, avg_F2, 'r-','LineWidth',1.5);
xlabel('迭代次数'); ylabel('第一前沿平均 F_2');
title('F_2 收敛曲线'); grid on;
saveas(gcf, 'Convergence_Curves.png');

%% ---------- 最优折衷解选择 (基于理想点距离) ----------
% 理想点: 各个目标的最小值
ideal_point = min(pareto_objectives);
% 归一化后计算距离
norm_pf = (pareto_objectives - ideal_point) ./ (max(pareto_objectives) - ideal_point + 1e-10);
distances = sqrt(sum(norm_pf.^2, 2));
[~, comp_idx] = min(distances);
compromise_solution = pareto_solutions(comp_idx, :);
compromise_obj = pareto_objectives(comp_idx, :);

fprintf('\n========== 最优折衷解 ==========\n');
fprintf('目标函数值: F1 = %.6f, F2 = %.6f\n', compromise_obj(1), compromise_obj(2));

% 提取 alpha, beta (度)
alpha_opt = compromise_solution(1:N_points);
beta_opt  = compromise_solution(N_points+1:end);

%% 图4: 最优折衷解的 alpha, beta 沿刀位点分布
figure('Color','white','Position',[100,100,800,600]);
subplot(2,1,1);
plot(1:N_points, alpha_opt, 'b-o','LineWidth',1.5,'MarkerSize',6);
hold on;
plot(1:N_points, lb_alpha, 'k--', 1:N_points, ub_alpha, 'k--');
xlabel('刀位点序号'); ylabel('\alpha (度)');
title('最优折衷解 \alpha 分布'); grid on;
legend('\alpha 最优值', '下界', '上界','Location','best');

subplot(2,1,2);
plot(1:N_points, beta_opt, 'r-o','LineWidth',1.5,'MarkerSize',6);
hold on;
plot(1:N_points, lb_beta, 'k--', 1:N_points, ub_beta, 'k--');
xlabel('刀位点序号'); ylabel('\beta (度)');
title('最优折衷解 \beta 分布'); grid on;
legend('\beta 最优值', '下界', '上界','Location','best');
saveas(gcf, 'Optimal_Alpha_Beta.png');

%% 图5: 最优折衷解相邻刀轴夹角 \gamma_i
[~, ~, gamma_opt] = calc_objectives(compromise_solution, theta_deg, RT, c0);
figure('Color','white','Position',[100,100,600,400]);
plot(2:N_points, gamma_opt, 's-','LineWidth',1.5,'MarkerSize',6,'Color',[0.8,0.2,0.2]);
xlabel('刀位点 i'); ylabel('\gamma_i (弧度)');
title('最优折衷解相邻刀轴向量夹角');
grid on;
saveas(gcf, 'Gamma_Angles.png');

%% 图6: Pareto 前沿上标注折衷解
figure('Color','white','Position',[100,100,800,600]);
plot(pareto_objectives(:,1), pareto_objectives(:,2), 'bo','MarkerSize',8);
hold on;
plot(compromise_obj(1), compromise_obj(2), 'rs','MarkerSize',12,'MarkerFaceColor','r');
xlabel('F_1'); ylabel('F_2');
title('Pareto 前沿与最优折衷解');
legend('Pareto 解','最优折衷解','Location','best');
grid on;
saveas(gcf, 'Pareto_with_Compromise.png');

%% ---------- 结果数据保存 ----------
save('Optimization_Results.mat', 'pareto_solutions', 'pareto_objectives', ...
    'compromise_solution', 'compromise_obj', 'alpha_opt', 'beta_opt', 'gamma_opt', ...
    'pareto_front_history', 'RT', 'c0');

fprintf('\n优化完成，所有图表与数据已保存。\n');

% =========================================================================
%                               NSWOA 主函数
% =========================================================================
function [final_pop, final_obj, pf_history] = NSWOA(obj_func, dim, lb, ub, ...
    pop_size, max_iter, n_obj, b_woa, theta_deg, RT, c0)
% 非支配排序鲸鱼优化算法 (NSWOA)
% 输入:
%   obj_func : 目标函数句柄 (实际通过 calc_objectives 实现)
%   dim      : 决策变量维数
%   lb, ub   : 变量边界 (1 x dim)
%   pop_size : 种群大小
%   max_iter : 最大迭代次数
%   n_obj    : 目标个数
%   b_woa    : 螺旋常数
%   theta_deg: 各点 theta 值 (度)
%   RT, c0   : 问题常数
% 输出:
%   final_pop : 最终种群决策变量
%   final_obj : 最终种群目标值
%   pf_history: 每代第一前沿目标值 (cell)

    % 初始化种群
    pop = repmat(lb, pop_size, 1) + repmat(ub - lb, pop_size, 1) .* rand(pop_size, dim);
    pop_obj = zeros(pop_size, n_obj);
    for i = 1:pop_size
        [f1, f2] = calc_objectives(pop(i,:), theta_deg, RT, c0);
        pop_obj(i,:) = [f1, f2];
    end

    % 非支配排序并选择初始种群
    [front_rank, crowd_dist] = NonDominatedSorting(pop_obj);
    pf_history = cell(max_iter, 1);

    % 主迭代循环
    for iter = 1:max_iter
        % 记录当前代的第一前沿
        first_front_idx = (front_rank == 1);
        pf_history{iter} = pop_obj(first_front_idx, :);

        % 选择领导者池 (第一前沿)
        leader_pool = pop(first_front_idx, :);
        leader_obj = pop_obj(first_front_idx, :);
        n_leaders = size(leader_pool, 1);

        if n_leaders == 0
            % 退化情况：随机选一个
            leader_pool = pop(randi(pop_size), :);
            n_leaders = 1;
        end

        % 子代种群预分配
        offspring = zeros(pop_size, dim);

        % WOA 参数
        a = 2 - iter * (2 / max_iter);          % 线性从 2 降至 0

        for i = 1:pop_size
            % 随机选择领导者
            leader = leader_pool(randi(n_leaders), :);

            % 更新参数
            r1 = rand(1, dim); r2 = rand(1, dim);
            A = 2 * a * r1 - a;
            C = 2 * r2;
            l = (a - 1) * rand(1, dim) + 1;    % l 在 [-1,1] 范围? 标准: l 在 [-1,1] 随机
            % 为保证 l ∈ [-1,1]，简单使用:
            l = 2 * rand(1, dim) - 1;

            p = rand();
            for j = 1:dim
                if p < 0.5
                    if abs(A(j)) < 1
                        % 收缩包围 (向领导者)
                        D = abs(C(j) * leader(j) - pop(i,j));
                        offspring(i,j) = leader(j) - A(j) * D;
                    else
                        % 随机搜索 (探索)
                        rand_idx = randi(pop_size);
                        X_rand = pop(rand_idx, j);
                        D = abs(C(j) * X_rand - pop(i,j));
                        offspring(i,j) = X_rand - A(j) * D;
                    end
                else
                    % 螺旋更新
                    dist = abs(leader(j) - pop(i,j));
                    offspring(i,j) = dist * exp(b_woa * l(j)) * cos(2 * pi * l(j)) + leader(j);
                end
                % 边界吸收
                offspring(i,j) = max(min(offspring(i,j), ub(j)), lb(j));
            end
        end

        % 计算子代目标函数
        offspring_obj = zeros(pop_size, n_obj);
        for i = 1:pop_size
            [f1, f2] = calc_objectives(offspring(i,:), theta_deg, RT, c0);
            offspring_obj(i,:) = [f1, f2];
        end

        % 合并父代与子代
        combined_pop = [pop; offspring];
        combined_obj = [pop_obj; offspring_obj];

        % 非支配排序并选择新一代
        [front_rank_comb, crowd_dist_comb] = NonDominatedSorting(combined_obj);

        % 精英选择: 逐层填充直到达到 pop_size
        new_pop = zeros(pop_size, dim);
        new_obj = zeros(pop_size, n_obj);
        count = 0;
        front = 1;
        while count + sum(front_rank_comb == front) <= pop_size
            idx = (front_rank_comb == front);
            n_add = sum(idx);
            new_pop(count+1:count+n_add, :) = combined_pop(idx, :);
            new_obj(count+1:count+n_add, :) = combined_obj(idx, :);
            count = count + n_add;
            front = front + 1;
        end
        % 若最后一层无法完全加入，按拥挤度排序选择部分个体
        if count < pop_size
            idx = find(front_rank_comb == front);
            % 按拥挤度距离降序排列
            [~, sort_idx] = sort(crowd_dist_comb(idx), 'descend');
            idx_sorted = idx(sort_idx);
            n_need = pop_size - count;
            new_pop(count+1:end, :) = combined_pop(idx_sorted(1:n_need), :);
            new_obj(count+1:end, :) = combined_obj(idx_sorted(1:n_need), :);
        end

        % 更新种群
        pop = new_pop;
        pop_obj = new_obj;
        [front_rank, crowd_dist] = NonDominatedSorting(pop_obj);

        % 打印进度
        if mod(iter, 20) == 0
            fprintf('迭代 %d/%d 完成, 第一前沿解数量: %d\n', ...
                iter, max_iter, sum(front_rank==1));
        end
    end

    final_pop = pop;
    final_obj = pop_obj;
end

% =========================================================================
%                          目标函数计算
% =========================================================================
function [F1, F2, gamma] = calc_objectives(x, theta_deg, RT, c0)
% 输入:
%   x        : 决策变量向量 [alpha(1..31), beta(1..31)] (度)
%   theta_deg: 各点 theta 值 (度)
%   RT, c0   : 问题常数
% 输出:
%   F1 : 目标函数 1 (负平均 I)
%   F2 : 目标函数 2 (平均 gamma^2)
%   gamma : 相邻刀轴夹角 (弧度) (用于后续分析)

    N = length(theta_deg);
    alpha_deg = x(1:N);
    beta_deg  = x(N+1:2*N);

    % 转换为弧度
    alpha = deg2rad(alpha_deg);
    beta  = deg2rad(beta_deg);
    theta = deg2rad(theta_deg);

    % 常数 R_G
    R_G = sqrt(RT^2 - (RT - c0)^2);

    % 目标函数 1
    B = (RT - c0) * cot(alpha);
    C = (RT - c0) * tan(beta - theta);
    I = pi * R_G^2 .* (R_G^2/2 + B.^2 + C.^2);
    F1 = -mean(I);

    % 目标函数 2: 计算 D_i 和 TA_i 向量
    D_val = sqrt(sin(alpha).^2 .* cos(beta).^2 + sin(beta).^2);
    TA = zeros(3, N);
    TA(1,:) = cos(alpha) .* sin(beta) ./ D_val;
    TA(2,:) = sin(alpha) .* cos(beta) ./ D_val;
    TA(3,:) = sin(alpha) .* sin(beta) ./ D_val;

    % 相邻夹角
    gamma = zeros(1, N-1);
    for i = 2:N
        dot_val = dot(TA(:,i-1), TA(:,i));
        dot_val = max(min(dot_val, 1), -1);  % 防止数值误差
        gamma(i-1) = acos(dot_val);
    end
    F2 = mean(gamma.^2);
end

% =========================================================================
%                    非支配排序与拥挤度计算
% =========================================================================
function [front_rank, crowd_dist] = NonDominatedSorting(obj_values)
% 对目标函数值 obj_values (n x n_obj) 进行非支配排序
% 返回每个个体的前沿等级 front_rank 和拥挤度距离 crowd_dist
% 基于 NSGA-II 方法

    [n, ~] = size(obj_values);
    front_rank = zeros(n, 1);
    crowd_dist = zeros(n, 1);
    
    % 支配关系矩阵
    dominated_count = zeros(n, 1);
    dominates_list = cell(n, 1);
    
    for i = 1:n
        for j = 1:n
            if i == j, continue; end
            if dominates(obj_values(i,:), obj_values(j,:))
                dominates_list{i} = [dominates_list{i}, j];
            elseif dominates(obj_values(j,:), obj_values(i,:))
                dominated_count(i) = dominated_count(i) + 1;
            end
        end
    end
    
    % 当前前沿
    current_front = find(dominated_count == 0);
    front = 1;
    while ~isempty(current_front)
        front_rank(current_front) = front;
        % 计算拥挤度距离
        crowd_dist(current_front) = crowding_distance(obj_values(current_front, :));
        
        next_front = [];
        for i = current_front(:)'
            for j = dominates_list{i}
                dominated_count(j) = dominated_count(j) - 1;
                if dominated_count(j) == 0
                    next_front = [next_front, j]; %#ok<AGROW>
                end
            end
        end
        front = front + 1;
        current_front = unique(next_front);
    end
end

function d = dominates(a, b)
% 判断 a 是否支配 b (最小化问题，a 所有目标 <= b 且至少一个严格 <)
    d = all(a <= b) && any(a < b);
end

function cd = crowding_distance(obj_sub)
% 计算子集 obj_sub 的拥挤度距离
    [m, n_obj] = size(obj_sub);
    cd = zeros(m, 1);
    if m <= 2
        cd(:) = Inf;
        return;
    end
    for j = 1:n_obj
        [~, idx] = sort(obj_sub(:,j));
        cd(idx(1)) = Inf;
        cd(idx(end)) = Inf;
        obj_range = obj_sub(idx(end), j) - obj_sub(idx(1), j);
        if obj_range == 0
            continue;
        end
        for k = 2:m-1
            cd(idx(k)) = cd(idx(k)) + ...
                (obj_sub(idx(k+1), j) - obj_sub(idx(k-1), j)) / obj_range;
        end
    end
end