% 定义 t 和 u 的范围
t = linspace(0, 10, 100); % t 的范围
u = linspace(0, 2, 100); % u 的范围

% 创建网格
[T, U] = meshgrid(t, u);

% 计算 f(t, u)
F = U .* (2 - U .* (sin(T)+ 2));

% 绘制三维曲面
figure;
surf(T, U, F);
xlabel('t');
ylabel('u');
zlabel('f(t, u)');
title('f(t, u) = u(1 - u sin(t))');
colorbar;
grid on;