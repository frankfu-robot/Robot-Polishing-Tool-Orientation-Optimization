clc;
clear;

%% ===================== 1. 常数设置 =====================
% 下面这些常数可以根据实际情况修改
RW    = 25;      % 工件半径或相关尺寸
RT    = 4;       % 刀具半径
C0    = 1;       % 常数 c0
l     = 60;      % 左侧边界或相关长度
LT    = 35;      % 刀具有效长度
dT    = 6;       % 刀具杆部半径或相关尺寸
dsafe = 1;     % 安全距离

%% ===================== 2. 刀位点坐标输入 =====================
% 每一行格式：[Pi, xi, yi, zi]
P = [
    1   20  -21.99999   0
    2   20  -21.8795   -2.2996
    3   20  -21.5192   -4.5741
    4   20  -20.9232   -6.7984
    5   20  -20.098    -8.9482
    6   20  -19.0526   -11
    7   20  -17.7984   -12.9313
    8   20  -16.3492   -14.7209
    9   20  -14.7209   -16.3492
    10  20  -12.9313   -17.7984
    11  20  -11        -19.0526
    12  20  -8.9482    -20.098
    13  20  -6.7984    -20.9232
    14  20  -4.5741    -21.5192
    15  20  -2.2996    -21.8795
    16  20   0         -22
    17  20   2.2996    -21.8795
    18  20   4.5741    -21.5192
    19  20   6.7984    -20.9232
    20  20   8.9482    -20.098
    21  20   11        -19.0526
    22  20   12.9313   -17.7984
    23  20   14.7209   -16.3492
    24  20   16.3492   -14.7209
    25  20   17.7984   -12.9313
    26  20   19.0526   -11
    27  20   20.098    -8.9482
    28  20   20.9232   -6.7984
    29  20   21.5192   -4.5741
    30  20   21.8795   -2.2996
    31  20   21.99999   0
];

Pi = P(:, 1);
x  = P(:, 2);
y  = P(:, 3);
z  = P(:, 4);

n = length(Pi);

%% ===================== 3. 结果变量初始化 =====================
theta = zeros(n, 1);   % 保存每个刀位点对应的 theta_i，单位为弧度

alpha_min_sph = zeros(n, 1);
alpha_max_sph = zeros(n, 1);
beta_min_sph  = zeros(n, 1);
beta_max_sph  = zeros(n, 1);

alpha_min_col = zeros(n, 1);
alpha_max_col = zeros(n, 1);
beta_min_col  = zeros(n, 1);
beta_max_col  = zeros(n, 1);

alpha_min_final = zeros(n, 1);
alpha_max_final = zeros(n, 1);
beta_min_final  = zeros(n, 1);
beta_max_final  = zeros(n, 1);

is_feasible = true(n, 1);

%% ===================== 4. sph 约束计算 =====================
% sph 约束中 alpha 的上下界与刀位点无关
ratio_sph = (RT - C0) / RT;

% 避免浮点误差导致 asin / acos 输入略微超出 [-1, 1]
ratio_sph = clampValue(ratio_sph, -1, 1);

alpha_sph_min_common = acos(ratio_sph);
alpha_sph_max_common = pi - acos(ratio_sph);

for i = 1:n

    xi = x(i);
    yi = y(i);
    zi = z(i);

    %% ---------- sph 约束下 alpha 区间 ----------
    alpha_min_sph(i) = alpha_sph_min_common;
    alpha_max_sph(i) = alpha_sph_max_common;

    %% ---------- sph 约束下 beta 区间 ----------
    theta_i = acos((-yi) / sqrt(yi^2 + zi^2));

    % 保存当前刀位点对应的 theta_i
    theta(i) = theta_i;

    beta_min_sph(i) = theta_i - asin(ratio_sph);
    beta_max_sph(i) = theta_i + asin(ratio_sph);

    %% ===================== 5. col 约束计算 =====================

    %% ---------- col 约束下 alpha 极小值 ----------
    PE_left = sqrt((l - xi)^2 + zi^2);

    if PE_left <= LT
        ratio = (RT + dsafe) / PE_left;
        alpha_min_col(i) = atan2(abs(zi), l - xi) + asinSafe(ratio, i, 'alpha_min_col');
    else
        ratio = (dT + dsafe) / PE_left;
        alpha_min_col(i) = atan2(abs(zi), l - xi) + asinSafe(ratio, i, 'alpha_min_col');
    end

    %% ---------- col 约束下 alpha 极大值 ----------
    PE_right = sqrt(xi^2 + zi^2);

    if PE_right <= LT
        ratio = (RT + dsafe) / PE_right;
        alpha_max_col(i) = pi - atan2(abs(zi), xi) - asinSafe(ratio, i, 'alpha_max_col');
    else
        ratio = (dT + dsafe) / PE_right;
        alpha_max_col(i) = pi - atan2(abs(zi), xi) - asinSafe(ratio, i, 'alpha_max_col');
    end

    %% ---------- col 约束下 beta 极小值 ----------
    PF_right = sqrt((RW - yi)^2 + zi^2);

    if PF_right <= (RT + dsafe)
        beta_min_col(i) = 0;
    elseif PF_right <= LT
        ratio = (RT + dsafe) / PF_right;
        beta_min_col(i) = atan2(abs(zi), RW - yi) + asinSafe(ratio, i, 'beta_min_col');
    else
        ratio = (dT + dsafe) / PF_right;
        beta_min_col(i) = atan2(abs(zi), RW - yi) + asinSafe(ratio, i, 'beta_min_col');
    end

    %% ---------- col 约束下 beta 极大值 ----------
    PF_left = sqrt((yi + RW)^2 + zi^2);

    if PF_left <= (RT + dsafe)
        beta_max_col(i) = pi;
    elseif PF_left <= LT
        ratio = (RT + dsafe) / PF_left;
        beta_max_col(i) = pi - atan2(abs(zi), yi + RW) - asinSafe(ratio, i, 'beta_max_col');
    else
        ratio = (dT + dsafe) / PF_left;
        beta_max_col(i) = pi - atan2(abs(zi), yi + RW) - asinSafe(ratio, i, 'beta_max_col');
    end

    %% ===================== 6. sph 与 col 区间求交集 =====================
    alpha_min_final(i) = max(alpha_min_sph(i), alpha_min_col(i));
    alpha_max_final(i) = min(alpha_max_sph(i), alpha_max_col(i));

    beta_min_final(i) = max(beta_min_sph(i), beta_min_col(i));
    beta_max_final(i) = min(beta_max_sph(i), beta_max_col(i));

    %% ---------- 判断当前刀位点是否存在可行姿态角区间 ----------
    if alpha_min_final(i) > alpha_max_final(i) || beta_min_final(i) > beta_max_final(i)
        is_feasible(i) = false;
    end

end

%% ===================== 7. 输出弧度制结果 =====================
Result_rad = table( ...
    Pi, x, y, z, theta, ...
    alpha_min_sph, alpha_max_sph, beta_min_sph, beta_max_sph, ...
    alpha_min_col, alpha_max_col, beta_min_col, beta_max_col, ...
    alpha_min_final, alpha_max_final, beta_min_final, beta_max_final, ...
    is_feasible);

disp('===================== 弧度制计算结果 =====================');
disp(Result_rad);

%% ===================== 8. 输出角度制结果 =====================
theta_deg = rad2deg(theta);

Result_deg = table( ...
    Pi, x, y, z, theta_deg, ...
    rad2deg(alpha_min_sph), rad2deg(alpha_max_sph), ...
    rad2deg(beta_min_sph),  rad2deg(beta_max_sph), ...
    rad2deg(alpha_min_col), rad2deg(alpha_max_col), ...
    rad2deg(beta_min_col),  rad2deg(beta_max_col), ...
    rad2deg(alpha_min_final), rad2deg(alpha_max_final), ...
    rad2deg(beta_min_final),  rad2deg(beta_max_final), ...
    is_feasible);

Result_deg.Properties.VariableNames = { ...
    'Pi', 'x', 'y', 'z', 'theta_deg', ...
    'alpha_min_sph_deg', 'alpha_max_sph_deg', ...
    'beta_min_sph_deg',  'beta_max_sph_deg', ...
    'alpha_min_col_deg', 'alpha_max_col_deg', ...
    'beta_min_col_deg',  'beta_max_col_deg', ...
    'alpha_min_final_deg', 'alpha_max_final_deg', ...
    'beta_min_final_deg',  'beta_max_final_deg', ...
    'is_feasible'};

disp('===================== 角度制计算结果 =====================');
disp(Result_deg);

%% ===================== 9. 保存结果到 Excel 文件 =====================
writetable(Result_rad, 'tool_orientation_range_rad.xlsx');
writetable(Result_deg, 'tool_orientation_range_deg.xlsx');

disp('===================== 文件保存完成 =====================');
disp('弧度制结果已保存为：tool_orientation_range_rad.xlsx');
disp('角度制结果已保存为：tool_orientation_range_deg.xlsx');

%% ===================== 10. 可选：单独输出每个刀位点的 theta_i =====================
Theta_Result = table(Pi, x, y, z, theta, theta_deg);

disp('===================== 每个刀位点对应的 theta_i =====================');
disp(Theta_Result);

writetable(Theta_Result, 'theta_result.xlsx');
disp('theta_i 结果已保存为：theta_result.xlsx');

%% ===================== 11. 局部函数 =====================

function value = clampValue(value, lower_bound, upper_bound)
    % 将数值限制在指定区间内
    % 主要用于避免 acos / asin 输入因为浮点误差略微超出 [-1, 1]

    value = max(lower_bound, min(upper_bound, value));
end

function angle = asinSafe(ratio, point_id, var_name)
    % 安全计算 asin
    % 如果 ratio 由于浮点误差略微超过 [-1, 1]，则进行截断
    % 如果 ratio 明显超界，则给出警告，并返回 NaN

    tol = 1e-10;

    if ratio > 1 && ratio <= 1 + tol
        ratio = 1;
    elseif ratio < -1 && ratio >= -1 - tol
        ratio = -1;
    elseif ratio > 1 || ratio < -1
        warning('第 %d 个刀位点的 %s 中 asin 输入 %.6f 超出 [-1, 1]，结果记为 NaN。', ...
            point_id, var_name, ratio);
        angle = NaN;
        return;
    end

    angle = asin(ratio);
end