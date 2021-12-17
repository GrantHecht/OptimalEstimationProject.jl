clear; close all; clc

%% Load data
load("./data/MATLABdata/ukf.mat", "em", "es", "s311", "s311mc", ...
    "s322", "s322mc", "s333", "s333mc", "s344", "s344mc", "s355", "s355mc", ...
    "s366", "s366mc", "s377", "s377mc", "ts", "xhat")

%% Position Errors
figure()
subplot(3,2,1)
plot(ts, es(:,1).*1000, 'r')
hold on
plot(ts, s311.*1000, 'Color', '#808080')
plot(ts, -s311.*1000, 'Color', '#808080')
plot(ts, s311mc.*1000, 'k')
plot(ts, -s311mc.*1000, 'k')
ylabel('$e_{r_x}$, m', 'Interpreter', 'latex')

set(gca, "fontname", "Times New Roman", "fontsize", 10)

xlim([0.0, 33.0])
ylim([-1000, 1000])
grid on

subplot(3,2,3)
plot(ts, es(:,1).*1000, 'r')
hold on
plot(ts, s322.*1000, 'Color', '#808080')
plot(ts, -s322.*1000, 'Color', '#808080')
plot(ts, s322mc.*1000, 'k')
plot(ts, -s322mc.*1000, 'k')
ylabel('$e_{r_y}$, m', 'Interpreter', 'latex')

set(gca, "fontname", "Times New Roman", "fontsize", 10)

xlim([0.0, 33.0])
ylim([-1500,1500])
grid on

subplot(3,2,5)
plot(ts, es(:,3).*1000, 'r')
hold on
plot(ts, s333.*1000, 'Color', '#808080')
plot(ts, -s333.*1000, 'Color', '#808080')
plot(ts, s333mc.*1000, 'k')
plot(ts, -s333mc.*1000, 'k')
xlabel('Time, days')
ylabel('$e_{r_z}$, m', 'Interpreter', 'latex')

set(gca, "fontname", "Times New Roman", "fontsize", 10)

xlim([0.0, 33.0])
ylim([-1000, 1000])
grid on

subplot(3,2,2)
plot(ts, es(:,4).*1000, 'r')
hold on
plot(ts, s344.*1000, 'Color', '#808080')
plot(ts, s344mc.*1000, 'k')
plot(ts, -s344.*1000, 'Color', '#808080')
plot(ts, -s344mc.*1000, 'k')
leg1 = legend('Estimation Error', 'UKF $3\sigma$', 'MC $3\sigma$','Interpreter', 'latex');
ylabel('$e_{v_x}$, m/s', 'Interpreter', 'latex')

set(gca, "fontname", "Times New Roman", "fontsize", 10)

ylim([-1, 1])
xlim([0.0, 33.0])
grid on

subplot(3,2,4)
plot(ts, es(:,5).*1000, 'r')
hold on
plot(ts, s355.*1000, 'Color', '#808080')
plot(ts, -s355.*1000, 'Color', '#808080')
plot(ts, s355mc.*1000, 'k')
plot(ts, -s355mc.*1000, 'k')
ylabel('$e_{v_y}$, m/s', 'Interpreter', 'latex')

set(gca, "fontname", "Times New Roman", "fontsize", 10)

ylim([-1.5, 1.5])
xlim([0.0, 33.0])
grid on

subplot(3,2,6)
plot(ts, es(:,6).*1000, 'r')
hold on
plot(ts, s366.*1000, 'Color', '#808080')
plot(ts, -s366.*1000, 'Color', '#808080')
plot(ts, s366mc.*1000, 'k')
plot(ts, -s366mc.*1000, 'k')
xlabel('Time, days')
ylabel('$e_{v_z}$, m/s', 'Interpreter', 'latex')

ylim([-1, 1])
xlim([0.0, 33.0])
grid on

set(gca, "fontname", "Times New Roman", "fontsize", 10)
set(leg1, "fontname", "Times New Roman", "fontsize", 8)
set(gcf, "PaperUnits","inches","PaperPosition",[0.25,0.25,6.5,8.0])
set(gcf, "PaperPositionMode","Manual")
print("./figures/UKFPosVelError.eps","-depsc","-r600")

%% Mass Error
figure()
plot(ts, es(:,7), 'r')
hold on
plot(ts, s377, 'Color', '#808080')
plot(ts, s377mc, 'k')
plot(ts, -s377, 'Color', '#808080')
plot(ts, -s377mc, 'k')
xlabel('Time, days')
ylabel('$e_m$, kg', 'Interpreter', 'latex')
leg2 = legend('Estimation Error', 'UKF $3\sigma$', 'MC $3\sigma$','Interpreter', 'latex');

grid on

set(gca, "fontname", "Times New Roman", "fontsize", 10)
set(leg2, "fontname", "Times New Roman", "fontsize", 8)
set(gcf, "PaperUnits","inches","PaperPosition",[0.25,0.25,6.5,2.7])
set(gcf, "PaperPositionMode","Manual")
print("./figures/UKFMassError.eps","-depsc","-r600")

%% rx Hist
figure()
subplot(3,2,1)
histogram(em(:,1).*1000)
hold on
xline(mean(em(:,1))*1000, '--k', 'LineWidth', 3)
xlabel('$e_{r_x}$, m', 'Interpreter', 'latex')
ylabel('# Occurances')
xlim([-10,10])
grid on

set(gca, "fontname", "Times New Roman", "fontsize", 10)

% ry Hist
subplot(3,2,3)
histogram(em(:,2).*1000)
hold on
xline(mean(em(:,2))*1000, '--k', 'LineWidth', 3)
xlabel('$e_{r_y}$, m', 'Interpreter', 'latex')
ylabel('# Occurances')
xlim([-15,15])
grid on

set(gca, "fontname", "Times New Roman", "fontsize", 10)

% rz Hist
subplot(3,2,5)
histogram(em(:,3).*1000)
hold on
xline(mean(em(:,3))*1000, '--k', 'LineWidth', 3)
xlabel('$e_{r_z}$, m', 'Interpreter', 'latex')
ylabel('# Occurances')
xlim([-10,10])
grid on

set(gca, "fontname", "Times New Roman", "fontsize", 10)

% vx Hist
subplot(3,2,2)
histogram(em(:,4).*1000)
hold on
xline(mean(em(:,4))*1000, '--k', 'LineWidth', 3)
xlabel('$e_{v_x}$, m/s', 'Interpreter', 'latex')
ylabel('# Occurances')
xlim([-0.01,0.01])
grid on

set(gca, "fontname", "Times New Roman", "fontsize", 10)

% vy Hist
subplot(3,2,4)
histogram(em(:,5).*1000)
hold on
xline(mean(em(:,5))*1000, '--k', 'LineWidth', 3)
xlabel('$e_{v_y}$, m/s', 'Interpreter', 'latex')
ylabel('# Occurances')
xlim([-0.02,0.02])
grid on

set(gca, "fontname", "Times New Roman", "fontsize", 10)

% vz Hist
subplot(3,2,6)
histogram(em(:,6).*1000)
hold on
xline(mean(em(:,6))*1000, '--k', 'LineWidth', 3)
xlabel('$e_{v_z}$, m/s', 'Interpreter', 'latex')
ylabel('# Occurances')
xlim([-0.01,0.01])
grid on

set(gca, "fontname", "Times New Roman", "fontsize", 10)
set(gcf, "PaperUnits","inches","PaperPosition",[0.25,0.25,6.5,8.0])
set(gcf, "PaperPositionMode","Manual")
print("./figures/UKFPosVelHist.eps","-depsc","-r600")

%% m Hist
figure()
histogram(em(:,7),60)
hold on
xline(mean(em(:,7)), '--k', 'LineWidth', 3)
xlabel('$e_m$, m/s', 'Interpreter', 'latex')
ylabel('# Occurances')
grid on

set(gca, "fontname", "Times New Roman", "fontsize", 10)
set(gcf, "PaperUnits","inches","PaperPosition",[0.25,0.25,3.25,2.7])
set(gcf, "PaperPositionMode","Manual")
print("./figures/UKFMassHist.eps","-depsc","-r600")