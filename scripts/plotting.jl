# Plotting
function plotEKF(ekf, xtrue, n)
    ts   = ekf.txp[1:n]
    xhat = ekf.xhats[1:n, :]
    xt   = xtrue[1:n, :]
    es   = ekf.es[1:n, :]
    σ311 = 3*sqrt.(ekf.Ps[1:n,1])
    σ322 = 3*sqrt.(ekf.Ps[1:n,2])
    σ333 = 3*sqrt.(ekf.Ps[1:n,3])
    σ344 = 3*sqrt.(ekf.Ps[1:n,4])
    σ355 = 3*sqrt.(ekf.Ps[1:n,5])
    σ366 = 3*sqrt.(ekf.Ps[1:n,6])
    σ377 = 3*sqrt.(ekf.Ps[1:n,7])

    mat"""
    ts      = $ts;
    ts      = (ts - $(ekf.scSim.t0))/86400;
    xhat    = $xhat;
    xt      = $xt;
    es      = $es;
    s311    = $σ311;
    s322    = $σ322;
    s333    = $σ333;
    s344    = $σ344;
    s355    = $σ355;
    s366    = $σ366;
    s377    = $σ377;

    figure()
    subplot(3,2,1)
    plot(ts, es(:,1), 'k')
    hold on
    plot(ts, s311, 'r')
    plot(ts, -s311, 'r')
    xlim([0, 33.2])
    ylim([-25, 25])
    grid on
    ylabel("\$r_x\$ Error, km", "Interpreter", "latex")

    subplot(3,2,3)
    plot(ts, es(:,1), 'k')
    hold on
    plot(ts, s322, 'r')
    plot(ts, -s322, 'r')
    xlim([0, 33.2])
    ylim([-100, 100])
    grid on
    ylabel("\$r_y\$ Error, km", "Interpreter", "latex")

    subplot(3,2,5)
    plot(ts, es(:,3), 'k')
    hold on
    plot(ts, s333, 'r')
    plot(ts, -s333, 'r')
    xlim([0, 33.2])
    ylim([-50, 50])
    xlabel("Time, days", "Interpreter", "latex")
    ylabel("\$r_z\$ Error, km", "Interpreter", "latex")
    grid on

    subplot(3,2,2)
    plot(ts, es(:,4)*1e3, 'k')
    hold on
    plot(ts, s344*1e3, 'r')
    plot(ts, -s344*1e3, 'r')
    xlim([0, 33.2])
    ylim([-6, 6])
    grid on
    ylabel("\$v_x\$ Error, m/s", "Interpreter", "latex")
    leg1 = legend('Estimation Error', 'EKF \$3\\sigma\$', 'Interpreter', 'latex');

    subplot(3,2,4)
    plot(ts, es(:,5)*1e3, 'k')
    hold on
    plot(ts, s355*1e3, 'r')
    plot(ts, -s355*1e3, 'r')
    xlim([0, 33.2])
    ylim([-10, 10])
    grid on
    ylabel("\$v_y\$ Error, m/s", "Interpreter", "latex")

    subplot(3,2,6)
    plot(ts, es(:,6)*1e3, 'k')
    hold on
    plot(ts, s366*1e3, 'r')
    plot(ts, -s366*1e3, 'r')
    xlim([0, 33.2])
    ylim([-6, 6])
    grid on
    xlabel("Time, days", "Interpreter", "latex")
    ylabel("\$v_z\$ Error, m/s", "Interpreter", "latex")

    set(gca, "fontname", "Times New Roman", "fontsize", 10)
    set(leg1, "fontname", "Times New Roman", "fontsize", 8)
    set(gcf, "PaperUnits", "inches", "PaperPosition", [0.25, 0.25, 8.0, 8.0])
    set(gcf, "PaperPositionMode", "Manual")
    print("./figures/EKFPosVelError.pdf", "-dpdf", "-r300")

    figure()
    plot(ts, es(:,7), 'k')
    hold on
    plot(ts, s377, 'r')
    plot(ts, -s377, 'r')
    xlim([0, 33.2])
    xlabel("Time, days", "Interpreter", "latex")
    ylabel("m Error, kg", "Interpreter", "latex")
    leg2 = legend('Estimation Error', 'EKF \$3\\sigma\$', 'Interpreter', 'latex');

    set(gca, "fontname", "Times New Roman", "fontsize", 10)
    set(leg2, "fontname", "Times New Roman", "fontsize", 8)
    set(gcf, "PaperUnits", "inches", "PaperPosition", [0.25, 0.25, 4.0, 3.0])
    set(gcf, "PaperPositionMode", "Manual")
    print("./figures/EKFMass.pdf", "-dpdf", "-r300")

    figure()
    plot3(xhat(:,1),xhat(:,2), xhat(:,3), 'r')
    hold on
    plot3(xt(:,1),xt(:,2),xt(:,3), 'b')
    grid on
    axis equal
    """
end

plotEKF(ekf, xtrue, length(ekf.txp))

# Create GPS plot
gpsTempSim = GPSSim(2170, 0, 2170, 1)
plotGPS(gpsTempSim)