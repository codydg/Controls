%% Run Simulation
%

%#ok<*UNRCH>
%clc, clear, close all;

GRAPHICAL_PLOT = false;

% Set up simulation parameters
M = 1000; % kg
m1 = 100; % kg
m2 = 100; % kg
l1 = 20; % m
l2 = 10; % m
g = 9.81; % m/s^2

% Create struct
params = struct;
params.M = M;
params.m1 = m1;
params.m2 = m2;
params.l1 = l1;
params.l2 = l2;
params.g = g;

stateLin = [0;0;0;0;pi/16;-pi/16];
stateNonLin = stateLin;
controlledStateLin = stateLin;
controlledStateNonLin = stateLin;

AF = [0,1,0,0,0,0;0,0,-g*m1/M,0,-g*m2/M,0;0,0,0,1,0,0;0,0,-g*(M+m1)/(M*l1),0,-g*m2/(M*l1),0;0,0,0,0,0,1;0,0,-g*m1/(M*l2),0,-g*(M+m2)/(M*l2),0];
BF = [0;1/M;0;1/(M*l1);0;1/(M*l2)];

Q = diag([1,1,1,1,1,1]);
R = 1;
K = lqr(AF,BF,Q,R);

% Set up Time data
step = 0.01; % Seconds
timesteps = 0:step:1000-step;
resultLin = zeros(numel(timesteps) + 1, 8);
resultNonLin = resultLin;
controlledResultLin = resultLin;
controlledResultNonLin = resultLin;
if GRAPHICAL_PLOT
    fig = figure; fig.Position = [17 100 1200 728];
    axLin = subplot(2,2,1); xlabel(axLin,'Lin');
    axNonLin = subplot(2,2,2); xlabel(axNonLin,'NonLin');
    axConLin = subplot(2,2,3); xlabel(axConLin,'Controlled Lin');
    axConNonLin = subplot(2,2,4); xlabel(axConNonLin,'Controlled NonLin');
end
for timeIndex = 1:numel(timesteps)
    FLin = 0;
    FNonLin = 0;
    FConLin = -K * controlledStateLin;
    FConNonLin = -K * controlledStateNonLin;
    
    resultLin(timeIndex,:) = [timesteps(timeIndex), stateLin.', FLin];
    resultNonLin(timeIndex,:) = [timesteps(timeIndex), stateNonLin.', FNonLin];
    controlledResultLin(timeIndex,:) = [timesteps(timeIndex), controlledStateLin.', FConLin];
    controlledResultNonLin(timeIndex,:) = [timesteps(timeIndex), controlledStateNonLin.', FConNonLin];
    
    stateLin = simulateLinearSystem(stateLin, FLin, step, params);
    stateNonLin = simulateNonLinearSystem(stateNonLin, FNonLin, step, params);
    controlledStateLin = simulateLinearSystem(controlledStateLin, FConLin, step, params);
    controlledStateNonLin = simulateNonLinearSystem(controlledStateNonLin, FConNonLin, step, params);
    
    if GRAPHICAL_PLOT
        plotState(axLin, stateLin(1), stateLin(3), stateLin(5), params.l1, params.l2);
        plotState(axNonLin, stateNonLin(1), stateNonLin(3), stateNonLin(5), params.l1, params.l2);
        plotState(axConLin, controlledStateLin(1), controlledStateLin(3), controlledStateLin(5), params.l1, params.l2);
        plotState(axConNonLin, controlledStateNonLin(1), controlledStateNonLin(3), controlledStateNonLin(5), params.l1, params.l2);
    end
end
resultLin(end,:) = [timesteps(end) + step, stateLin.', nan];
resultNonLin(end,:) = [timesteps(end) + step, stateNonLin.', nan];
controlledResultLin(end,:) = [timesteps(end) + step, controlledStateLin.', nan];
controlledResultNonLin(end,:) = [timesteps(end) + step, controlledStateNonLin.', nan];
figure('units','normalized','outerposition',[0 0 1 1]);
subplot 321;plot(resultLin(:,1),resultLin(:,2),'r',resultNonLin(:,1),resultNonLin(:,2),'b-.'); legend('Linear','Non-Linear'); ylabel('X'); xlabel('Time'); title('Uncontrolled X');
subplot 323;plot(resultLin(:,1),resultLin(:,4),'r',resultNonLin(:,1),resultNonLin(:,4),'b-.'); legend('Linear','Non-Linear'); ylabel('Theta 1'); xlabel('Time'); title('Uncontrolled Theta 1');
subplot 325;plot(resultLin(:,1),resultLin(:,6),'r',resultNonLin(:,1),resultNonLin(:,6),'b-.'); legend('Linear','Non-Linear'); ylabel('Theta 2'); xlabel('Time'); title('Uncontrolled Theta 2');
subplot 322;plot(controlledResultLin(:,1),controlledResultLin(:,2),'r',controlledResultNonLin(:,1),controlledResultNonLin(:,2),'b-.'); legend('Linear','Non-Linear'); ylabel('X'); xlabel('Time'); title('Controlled X');
subplot 324;plot(controlledResultLin(:,1),controlledResultLin(:,4),'r',controlledResultNonLin(:,1),controlledResultNonLin(:,4),'b-.'); legend('Linear','Non-Linear'); ylabel('Theta 1'); xlabel('Time'); title('Controlled Theta 1');
subplot 326;plot(controlledResultLin(:,1),controlledResultLin(:,6),'r',controlledResultNonLin(:,1),controlledResultNonLin(:,6),'b-.'); legend('Linear','Non-Linear'); ylabel('Theta 2'); xlabel('Time'); title('Controlled Theta 2');
function plotState(ax, x, theta1, theta2, l1, l2)
    axes(ax); hold off;
    
    % Plot Box
    boxWidth = 5; boxHeight = 5;
    plot(x + 0.5*[-boxWidth,boxWidth,boxWidth,-boxWidth,-boxWidth], ...
        0.5*[boxHeight,boxHeight,-boxHeight,-boxHeight,boxHeight],'k');
    hold on;
    
    % Plot Circles
    massRadius = 2.5;
    center1 = [x - l1*sin(theta1),-l1*cos(theta1)];
    viscircles(center1, massRadius, 'Color', 'b');
    center2 = [x - l2*sin(theta2),-l2*cos(theta2)];
    viscircles(center2, massRadius, 'Color', 'r');
    
    % Plot Ropes
    plot([center1(1), x, center2(1)], [center1(2), 0, center2(2)], 'k');
    
    % Adjust Axis
    axis([-60,60,-25,5]);
    drawnow;
end