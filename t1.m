% EXPERIMENT 1 – Asymmetric Vortex Swirl Thruster (plot stays open!)
clear all; close all; clc;

L  = 0.10;  R  = 0.020;  N  = 120;
I0 = 4.0;   f  = 6800;   omega = 2*pi*f;
eta_aether = 1.8e-18;

pitch_angles = [30, -30, 0, 90];
labels = {"+30° right-handed", "-30° left-handed", ...
          "0° symmetric cylinder", "90° almost straight"};
colors = {'r', 'b', 'k', 'm'};

n_seg = 600;
theta = linspace(0, N*2*pi, n_seg+1); theta(end)=[];
z     = linspace(0, L,      n_seg+1); z(end)=[];

res  = 30;  side = 0.65;
[xg,yg,zg] = meshgrid(linspace(-side/2, side/2, res));
dV   = (side/res)^3;

t = linspace(0, 5/f, 300);
thrust_all = zeros(length(t), 4);

fprintf('Running 4 cases (plot will stay open when finished)...\n');

for cas = 1:4
    fprintf('  Case %d/4: %s ... ', cas, labels{cas});
    pitch = pitch_angles(cas)*pi/180;
    drift = tan(pitch)*R*(z/L - 0.5);
    x = R*cos(theta);
    y = R*sin(theta) + drift;
    wire = [x(:) y(:) z(:)];

    thrust = zeros(size(t));
    for k = 1:length(t)
        if mod(k,60)==0, fprintf('.'), end
        I = I0*sin(omega*t(k));

        A = zeros(res,res,res,3);
        for i = 1:n_seg-1
            dl = wire(i+1,:) - wire(i,:);
            rm = wire(i,:) + 0.5*dl;
            rx = xg-rm(1); ry = yg-rm(2); rz = zg-rm(3);
            r3 = (rx.^2 + ry.^2 + rz.^2 + 1e-20).^(1.5);

            cross_x = dl(2)*rz - dl(3)*ry;
            cross_y = dl(3)*rx - dl(1)*rz;
            cross_z = dl(1)*ry - dl(2)*rx;

            contrib = 1e-7 * I ./ r3;
            A(:,:,:,1) += contrib .* cross_x;
            A(:,:,:,2) += contrib .* cross_y;
            A(:,:,:,3) += contrib .* cross_z;
        end

        F = -eta_aether * sum(A,[1 2 3]) * dV;
        thrust(k) = norm(F)*1e6;        % μN
    end
    thrust_all(:,cas) = thrust;
    fprintf(' done → %.1f μN\n', mean(thrust(50:end)));
end

% ====================== PLOT THAT STAYS OPEN ======================
figure('Name','Asymmetric Vortex Thruster – Fluid Hypothesis','NumberTitle','off');
hold on; grid on; box on;

for cas = 1:4
    plot(t*1e3, thrust_all(:,cas), 'LineWidth', 3, 'Color', colors{cas}, ...
         'DisplayName', labels{cas});
end

xlabel('Time (ms)', 'FontSize', 14);
ylabel('Predicted thrust (μN)', 'FontSize', 14);
title('Fluid Hypothesis – Asymmetric Vortex Swirl Thruster','FontSize', 16);
legend('Location','northwest','FontSize', 13);
xlim([0 t(end)*1e3]);
ylim([0 2400]);

% This single line is the magic → keeps the window open until you close it
disp('Simulation finished. Close the figure window when you are done looking.');
uiwait(gcf);       % ← pauses execution until you manually close the plot
% =================================================================

fprintf('\nWindow closed – goodbye!\n');
