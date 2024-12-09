close all; clc; clear;

load 'LRV_wheel_FRF.mat'

N = size(frf, 2);

[FRFopt, freqsOpt, xSol] = FRF_fit(freq, frf, 4, 200, 'complexDiff', true);

% Compare the experimental and identified mode shapes
theta_r = (0:11)*15;
thetaInterp = linspace(0, 11*15, 1000);
figure;
ax = polaraxes;
for i = 1:4
    A_r = xSol(3:N+2, i)/max(abs(xSol(3:N+2, i)));
    AInterp = interp1(theta_r, A_r, thetaInterp, 'spline');
    Aplot = polarplot(deg2rad(theta_r'), A_r/3+1, 'o', 'DisplayName', ['Mode ', num2str(i), ' - ', num2str(xSol(1, i)/(2*pi),4), ' Hz'], 'LineWidth', 1.5);
    hold on;
    interPlot = polarplot(deg2rad(thetaInterp), AInterp/3+1, 'DisplayName', ['Interpolation Mode ', num2str(i)], 'LineWidth', 1.5);
    interPlot.Color = Aplot.Color;
end
title('Mode shapes');
polarplot(deg2rad(0:1:360), ones(1, 361), 'k--', 'DisplayName', 'Original shape', 'LineWidth', 1.5);
legend();
ax.ThetaZeroLocation = 'bottom';
grid on;