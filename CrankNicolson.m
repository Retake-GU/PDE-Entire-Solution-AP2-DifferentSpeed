clear; clc;

%% 参数设置
L = 40;              % 空间范围 [-L, L]
dz = 0.1;            % 空间步长
z = -L:dz:L;         % 移动坐标系下的空间变量
Nz = length(z);

T_total = 100;       % 模拟总时间（反向时间从0到-100）
dt = 0.001;          % 时间步长
Nt = round(T_total/dt);
t_start = 0;         % 起始时间为0

c = 3;               % traveling wave 速度
T_period = 2*pi;     % 时间周期

% 非线性项
% a1 = @(t) 1.5 + sin(2*pi*t/T_period);
% f1 = @(t, v) v .* (a1(t) - v);
% P = @(t) (exp(1.5 * t - cos(t))) ./ (integral(@(s) exp(1.5 * s - cos(s)), 0, t) + 1);

% f1 = @(t, v) v .* (1 - a1(t) * v);
% P = @(t) 2 ./ (4 + sin(t) - cos(t) + exp(-t));

a1 = @(t) 2 + cos(2*pi*t/T_period);
f1 = @(t, v) v .* (3 - a1(t) * v);
P = @(t) 2 ./ (12 + sin(t) + cos(t) + 2 * exp(-t));
%% 初始条件
v = 1 ./ (1 + exp(z));    % 初始 traveling wave
v_next = zeros(size(v));

% 存储用于重建 u(t,x)（改为列向量）
z_all = [];
t_all = [];
v_all = [];

plot_interval = floor(Nt / 20);

%% Crank-Nicolson 方法参数
alpha = dt / (2 * dz^2);  % 隐式扩散项系数

%% 时间推进循环（反向时间）
for step = 1:Nt
    t = t_start - step * dt;  % 当前时间为负向推进
    
    % --- 显式计算对流项 (dvz) 和源项 (f1) ---
    dvz = zeros(Nz, 1);
    f1_val = zeros(Nz, 1);
    for i = 2:Nz-1
        dvz(i) = (v(i+1) - v(i-1)) / (2*dz);  % 二阶中心差分
        f1_val(i) = f1(t, v(i));              % 非线性源项
    end
    
    % --- 构建隐式线性系统 A * v_next = b ---
    % 主对角线 (i)
    main_diag = (1 + 2*alpha) * ones(Nz, 1);
    % 次对角线 (i-1 和 i+1)
    off_diag = -alpha * ones(Nz, 1);
    
    % 构造稀疏矩阵 A
    A = spdiags([off_diag, main_diag, off_diag], [-1, 0, 1], Nz, Nz);
    
    % 强制边界条件
    A(1, 1:2) = [1, 0];      % 左边界: v(1) = P(t)
    A(Nz, Nz-1:Nz) = [0, 1]; % 右边界: v(end) = 0
    
    % 构建右端向量 b
    b = zeros(Nz, 1);
    b(1) = P(t);             % 左边界值
    b(Nz) = 0;               % 右边界值
    
    % 显式部分填充 b
    for i = 2:Nz-1
        explicit_part = v(i) + dt/2 * ( (v(i+1) - 2*v(i) + v(i-1))/dz^2 + c*dvz(i) + f1_val(i) );
        b(i) = explicit_part;
    end
    
    % --- 求解线性系统 ---
    v_next = A \ b;
    
    % --- 更新解 ---
    v = v_next;
    
    % --- 记录数据（改为列向量拼接）---
    if mod(step, plot_interval) == 0
        z_all = [z_all; z(:)];          % 将 z 转为列向量
        t_all = [t_all; t * ones(Nz, 1)]; % 列向量
        v_all = [v_all; v(:)];          % 列向量
    end
end

%% 重建和绘图
% 生成网格点
z_vals = linspace(min(z_all), max(z_all), 200);
t_vals = linspace(min(t_all), max(t_all), 200);
[Z, T] = meshgrid(z_vals, t_vals);

% 插值（确保输入是列向量）
V = griddata(z_all, t_all, v_all, Z, T, 'cubic');

% 绘图
figure;
surf(Z, T, V, 'EdgeColor', 'none');
xlabel('z'); ylabel('t'); zlabel('phi_{1}^{c}(t,z)');
title('Crank-Nicolson: phi_1(t,z)');
colorbar;
view(45, 30);
grid on;

%% 截取 t = 0.2 时刻的二维图像
t_target = 0.2;
[~, idx] = min(abs(t_vals - t_target));
v_slice = V(idx, :);

figure;
plot(z_vals, v_slice, 'b-', 'LineWidth', 2);
xlabel('Transformed x');
ylabel(['u(t = ', num2str(t_vals(idx), '%.2f'), ', x)']);
title(['Crank-Nicolson Slice at t = ', num2str(t_vals(idx), '%.2f')]);
grid on;