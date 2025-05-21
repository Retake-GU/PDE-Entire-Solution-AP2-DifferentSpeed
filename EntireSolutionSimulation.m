clc; clear;

% 空间和时间参数
L = 50; Nx = 500;
x = linspace(0, L, Nx)';
dx = x(2) - x(1);

T_start = -50;
T_end = 50;
Nt = 1000;
dt = (T_end - T_start) / Nt;
t_vec = linspace(T_start, T_end, Nt);

% 初始条件
% u = exp(-((x - L/2).^2));
u = ones(size(x));
U = zeros(Nx, Nt);
U(:,1) = u;

% 定义 nu(t), f1(t,u), f2(t,u)
nu = @(t) exp(t) ./ (exp(t) + exp(-t));

T_period = 2*pi;
% a1 = @(t) 1 + sin(2*pi*t / T_period);
% a2 = @(t) 1 + cos(2*pi*t / T_period);
a1 = @(t)sin(2*pi*t / T_period);
a2 = @(t)sin(2*pi*t / T_period)/2;
% b1 = 1;
% b2 = 1;
b1 = @(t) cos(2*pi*t / T_period);
b2 = @(t) cos(2*pi*t / T_period)/2;

f1 = @(t, u) a1(t) .* u - b1(t) .* u.^2;
f2 = @(t, u) a2(t) .* u - b2(t) .* u.^2;
f = @(t, u) (1 - nu(t)) .* f1(t, u) + nu(t) .* f2(t, u);

% Crank–Nicolson 离散矩阵
alpha = dt / (2 * dx^2);
e = ones(Nx, 1);
A = spdiags([-alpha*e, (1 + 2*alpha)*e, -alpha*e], -1:1, Nx, Nx);
B = spdiags([alpha*e, (1 - 2*alpha)*e, alpha*e], -1:1, Nx, Nx);

% Dirichlet 边界条件处理
A(1,:) = 0; A(end,:) = 0; A(1,1) = 1; A(end,end) = 1;
B(1,:) = 0; B(end,:) = 0; B(1,1) = 1; B(end,end) = 1;

% 时间推进
for n = 1:Nt-1
    t = T_start + (n-1)*dt;
    F = f(t, u);
    b = B * u + dt * F;
    b(1) = 0; b(end) = 0;
    u = A \ b;
    U(:, n+1) = u;
end

% ===== 三维图：surf(t,x,u) =====
[T_grid, X_grid] = meshgrid(t_vec, x);

figure;
surf(T_grid, X_grid, U, 'EdgeColor', 'none');
xlabel('Time t');
ylabel('Space x');
zlabel('u(t,x)');
title('Solution u(t,x) of the PDE over t ∈ [-10,10]');
colorbar;
view(45, 30);
