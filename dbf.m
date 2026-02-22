% DBF Performance Analysis: Lift and Drag Buildup
clear; clc; close all;

%% --- INPUT DATA ---
W = 15.53;          % Mission 2 Weight (lb) - Baseline for plots
S_in2 = 590;        % Wing Area (sq in)
S = S_in2 / 144;    % Area in sq ft
AR = 5.9;           % Aspect Ratio
e = 0.80;           % Oswald Efficiency (Composite)
Cd0 = 0.032;        % Parasitic Drag Estimate
rho = 0.002377;     % Sea level density
V = 40:1:170;       % Speed range (ft/s)

%% --- CALCULATIONS ---
K = 1 / (pi * AR * e);
Cl = (2 * W) ./ (rho * S * V.^2);
Cd = Cd0 + K * Cl.^2;

% Drag Forces (lbf)
D_parasitic = 0.5 * rho * V.^2 * S * Cd0;
D_induced = 0.5 * rho * V.^2 * S .* (K * Cl.^2);
D_total = D_parasitic + D_induced;

%% --- PLOTTING ---
figure('Color', 'w', 'Position', [100, 100, 1100, 480]);

% --- LEFT PLOT: COEFFICIENTS ---
subplot(1,2,1);
yyaxis left
plot(V, Cl, 'LineWidth', 2, 'Color', [0 0.447 0.741]);
ylabel('Coefficient of Lift (C_L)'); ylim([0 1.5]);
yyaxis right
plot(V, Cd, 'LineWidth', 2, 'Color', [0.85 0.325 0.098]);
ylabel('Coefficient of Drag (C_D)'); ylim([0 0.15]);
grid on; xlabel('Speed (ft/s)');
title('Lift and Drag Coefficient');
legend('Required C_L', 'Coefficient of Drag', 'Location', 'northeast');

% --- RIGHT PLOT: FORCES ---
subplot(1,2,2);
plot(V, D_parasitic, '--', 'LineWidth', 1.5, 'Color', [0 0.447 0.741]); hold on;
plot(V, D_induced, '--', 'LineWidth', 1.5, 'Color', [0.929 0.694 0.125]);
plot(V, D_total, 'Color', [0.494 0.184 0.556], 'LineWidth', 2);
grid on; xlabel('Speed (ft/s)'); ylabel('Force (lbf)');
title('Lift and Drag Force Buildup');
legend('Parasitic Drag', 'Induced Drag', 'Total Drag', 'Location', 'northeast');

% Optional: Add marker for M2 Cruise Speed
plot(92, 1.714, 'ko', 'MarkerFaceColor', 'r'); % M2 Operating Point
text(94, 1.8, 'M2 Cruise', 'FontSize', 9);