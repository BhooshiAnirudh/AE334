clc; clear;

% Given parameters
a = 0.5;  % length (m)
b = 0.5;  % breadth (m)
t = 4e-3; % Thickness (m)
q0 = 1000; % Transverse load (N/m^2)
E = 200e9; % Young's Modulus (Pa)
neu = 0.3; % Poisson's ratio

% Flexural rigidity
D = E * (t^3) / (12 * (1 - neu^2));

% Symbolic variables for x and y
syms x y

% Hermite cubic basis function
lex = a; % Characteristic length in x
ley = b; % Characteristic length in y

phiX = [(x - a)^2 / lex^2 * (1 + 2*(x - 0)/lex);
         (x - a)^2 / lex^2 * (x - 0);
         (x - 0)^2 / lex^2 * (1 - 2*(x - a)/lex);
         (x - 0)^2 / lex^2 * (x - a)];

phiY = [(y - b)^2 / ley^2 * (1 + 2*(y - 0)/ley);
         (y - b)^2 / ley^2 * (y - 0);
         (y - 0)^2 / ley^2 * (1 - 2*(y - b)/ley);
         (y - 0)^2 / ley^2 * (y - b)];

% Basis functions
phi = cell(16,1);
for i = 1:4
    for j = 1:4
        index = 4*(i - 1) + j;
        phi{index} = phiX(i) * phiY(j);
    end
end

% Second derivatives
phiXX = diff(phi, x, 2);
phiYY = diff(phi, y, 2);
phiXY = diff(diff(phi, x), y);

% Stiffness matrix
K = sym(zeros(16,16));
for i = 1:16
    for j = 1:16
        K(i,j) = D * int(int( ...
            phiXX(j) * phiXX(i) + phiYY(j) * phiYY(i) + ...
            neu * (phiXX(j) * phiYY(i) + phiYY(j) * phiXX(i)) + ...
            2*(1-neu) * phiXY(j) * phiXY(i), ...
            x, 0, a), y, 0, b);
    end
end

% Force vector
F = sym(zeros(16,1));
for i = 1:16
    F(i) = int(int(q0 * phi{i}, x, 0, a), y, 0, b);
end

nonzeroAlphas = [6, 8, 14, 16];

% Boundary conditions
Kmod = K;
Fmod = F;
alpha = sym('alpha', [16,1]);

for k = 1:16
    if ~ismember(k, nonzeroAlphas)
        alpha(k) = 0;
        Fmod = Fmod - Kmod(:,k) * alpha(k);
        Kmod(:, k) = 0; 
        Kmod(k, :) = 0; 
        Kmod(k, k) = 1;
        Fmod(k) = alpha(k);
    end
end

% Convert to numerical values
KmodNum = double(Kmod); 
FmodNum = double(Fmod); 

% Solve for unknown alpha values
solNum = KmodNum \ FmodNum;

% Display results
disp(['Rank of Kmod: ', num2str(rank(KmodNum))]);
disp('Stiffness Matrix Kmod:');
disp(vpa(KmodNum,4)); 

disp('Force Vector F:');
disp(vpa(FmodNum,4));

disp('Alpha values:');
for k = 1:16
    fprintf('alpha(%d) = %.10f\n', k, solNum(k));
end

% Deflection function w(x,y)
wXY = 0;
for k = 1:16
    wXY = wXY + solNum(k) * phi{k};
end

% Strain energy
strainEnergy = 0.5 * sum(solNum .* (K * solNum));
disp(['Strain Energy: ', num2str(double(strainEnergy))]);

wFunc = matlabFunction(wXY, 'Vars', [x, y]);

% Grid for contour plot
[X, Y] = meshgrid(linspace(0, a, 50), linspace(0, b, 50));
W = wFunc(X, Y);

% Contour plot
figure;
contourf(X, Y, W, 20, 'LineColor', 'none');
colorbar;
xlabel('x (m)');
ylabel('y (m)');
title('w(x, y)');