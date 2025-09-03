clc; clear;

% Given parameters
a = 0.5;  % Length (m)
b = 0.5;  % Breadth (m)
t = 4e-3; % Thickness (m)
q0 = 1000; % Transverse load (N/m^2)
E = 200e9; % Young's Modulus (Pa)
neu = 0.3; % Poisson's ratio

% Flexural rigidity
D = E * (t^3) / (12 * (1 - neu^2));

% Symbolic variables for x and y
syms x y

% Hermite cubic basis functions
lEx = a; % Characteristic length in x
lEy = b; % Characteristic length in y

phiX = [((x - a)^2 / lEx^2) * (1 + (2/lEx) * (x - 0));
         ((x - a)^2 / lEx^2) * ((x - 0));
         ((x - 0)^2 / lEx^2) * (1 - (2/lEx) * (x - a));
         (1 / lEx^2) * ((x - 0)^2) * (x - a)];

phiY = [((y - b)^2 / lEy^2) * (1 + (2/lEy) * (y - 0));
         ((y - b)^2 / lEy^2) * ((y - 0));
         ((y - 0)^2 / lEy^2) * (1 - (2/lEy) * (y - b));
         (1 / lEy^2) * ((y - 0)^2) * (y - b)];

% Basis functions
phiC = sym(zeros(16,1));
for i = 1:4
    for j = 1:4
        I = 4*(i-1) + j; % Global index
        phiC(I) = phiX(i) * phiY(j);
    end
end

% Second derivatives
phiC_xx = diff(phiC, x, 2);
phiC_yy = diff(phiC, y, 2);
phiC_xy = diff(diff(phiC, x), y);

% Stiffness matrix
K = sym(zeros(16,16));
for i = 1:16
    for j = 1:16
        K(i,j) = D * int(int( ...
            phiC_xx(j) * phiC_xx(i) + phiC_yy(j) * phiC_yy(i) + ...
            neu * (phiC_xx(j) * phiC_yy(i) + phiC_yy(j) * phiC_xx(i)) + ...
            2*(1-neu) * phiC_xy(j) * phiC_xy(i), ...
            x, 0, a), y, 0, b);
    end
end

% Force vector
F = sym(zeros(16,1));
for i = 1:16
    F(i) = int(int(q0 * phiC(i), x, 0, a), y, 0, b);
end

nonzeroAlphas = [3, 4, 7, 8, 11, 12, 15, 16]; % These alphas remain unknown

% Boundary conditions
K_mod = K;
F_mod = F;
alpha = sym('alpha', [16,1]); % Symbolic alpha vector

% Apply boundary conditions: Set alpha(k) = 0 for specified indices
for k = 1:16
    if ~ismember(k, nonzeroAlphas)  % If alpha(k) must be set to zero
        alpha(k) = 0;
        F_mod = F_mod - K_mod(:,k) * alpha(k); % Update force vector
        K_mod(:, k) = 0; % Zero out column k
        K_mod(k, :) = 0; % Zero out row k
        K_mod(k, k) = 1; % Identity-like modification
        F_mod(k) = alpha(k); % Ensure right-hand side is zero
    end
end

% Convert to numerical values
K_mod_num = double(K_mod); 
F_mod_num = double(F_mod); 

% Solve for unknown alpha values
sol_num = K_mod_num \ F_mod_num;

% Display results
disp(['Rank of K_mod: ', num2str(rank(K_mod_num))]);
disp('Stiffness Matrix K_mod:');
disp(vpa(K_mod_num,4)); % Display with 4 decimal places

disp('Force Vector F:');
disp(vpa(F_mod_num,4));

disp('Alpha values:');
for k = 1:16
    fprintf('alpha(%d) = %.10f\n', k, sol_num(k));
end

% Deflection function w(x,y)
w_xy = sum(sol_num .* phiC);

w_func = matlabFunction(w_xy, 'Vars', [x, y]);

% Grid for contour plot
[X, Y] = meshgrid(linspace(0, a, 50), linspace(0, b, 50));
W = w_func(X, Y);

disp('Equation for deflection w(x,y):');
disp(vpa(w_xy, 4));

% Contour plot
figure;
contourf(X, Y, W, 20, 'LineColor', 'none');
colorbar;
xlabel('x (m)');
ylabel('y (m)');
title('w(x, y)');

% Strain energy
strainEnergy = 0.5 * int(int( ...
    D * (phiC_xx.' * K_mod * phiC_xx + phiC_yy.' * K_mod * phiC_yy + ...
    2*(1-neu) * phiC_xy.' * K_mod * phiC_xy), x, 0, a), y, 0, b);

disp('Strain Energy:');
disp(vpa(strainEnergy, 6));

% Basis functions
for I = 1:16
    fprintf('Phi%d(x,y) = %s\n', I, char(phiC(I)));
end