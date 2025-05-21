% 模拟 Fisher-KPP 类型的时间周期 traveling wave 解
clear; clc;

%% 参数设置
L = 20;              % 空间范围 [-L, L]
dz = 0.1;            % 空间步长
z = -L:dz:L;         % 空间网格
Nz = length(z);

T_total = 2;         % 模拟总时间
dt = 0.001;          % 时间步长
Nt = round(T_total/dt);

c = 2;               % 行波速度
T_period = 1;        % 时间周期

%% 初始条件（sigmoid 形状）
v = 1 ./ (1 + exp(z));    % 近似 traveling wave
v_next = zeros(size(v));

%% 函数定义
% b1 = @(t) 1 + sin(2*pi*t/T_period);
% f = @(t, v) v .* (b1(t) - v);
% P = @(t) b1(t);  % 左边界值 P_i(t)
% 改为 f2(t,v)
b2 = @(t) 1 + cos(2*pi*t/T_period);
f = @(t, v) v .* (b2(t) - v);
P = @(t) b2(t);
%% 可视化设置
figure;
plot_interval = floor(Nt/10);
plot_count = 0;

%% 时间推进
for step = 1:Nt
    t = step * dt;

    % 内部点更新（中心差分 + 显式 Euler）
    for i = 2:Nz-1
        d2vz = (v(i+1) - 2*v(i) + v(i-1)) / dz^2;
        dvz  = (v(i+1) - v(i-1)) / (2*dz);
        v_next(i) = v(i) + dt * (d2vz + c * dvz + f(t, v(i)));
    end

    % 边界条件
    v_next(1) = P(t);  % 左端趋近 P_i(t)
    v_next(end) = 0;   % 右端为 0

    % 更新解
    v = v_next;

    % 绘图
    if mod(step, plot_interval) == 0
        plot(z, v, 'DisplayName', sprintf('t = %.2f', t)); hold on;
        plot_count = plot_count + 1;
    end
end

%% 显示图像
xlabel('z'); ylabel('v(t,z)');
title('Time-periodic Traveling Wave Simulation (MATLAB)');
legend show;
grid on;
