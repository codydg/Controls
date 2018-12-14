%% Run Simulation
%

%#ok<*UNRCH>
clc, clear, close all;

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

xLin = 0;%input('Initial X: ');
theta1Lin = 0;%input('Initial Theta 1: ');
theta2Lin = 0;%input('Initial Theta 2: ');
dxLin = 0;%input('Initial d/dt X: ');
dtheta1Lin = pi/4;%input('Initial d/dt Theta 1: ');
dtheta2Lin = -pi/4;%input('Initial d/dt Theta 2: ');

xNonLin = xLin;
theta1NonLin = theta1Lin;
theta2NonLin = theta2Lin;
dxNonLin = dxLin;
dtheta1NonLin = dtheta1Lin;
dtheta2NonLin = dtheta2Lin;

AF = [0,1,0,0,0,0;0,0,-g*m1/M,0,-g*m2/M,0;0,0,0,1,0,0;0,0,-g*(M+m1)/(M*l1),0,-g*m2/(M*l1),0;0,0,0,0,0,1;0,0,-g*m1/(M*l2),0,-g*(M+m2)/(M*l2),0];
BF = [0;1/M;0;1/(M*l1);0;1/(M*l2)];

Q = diag([1,1,1,1,1,1]);
R = 1;
K = lqr(AF,BF,Q,R);

% Set up Time data
step = 0.3; % Seconds
timesteps = 0:step:100-step;
resultLin = zeros(numel(timesteps) + 1, 8);
resultNonLin = resultLin;
if GRAPHICAL_PLOT
    fig = figure; fig.Position = [17 100 1200 728];
    axLin = subplot(2,1,1);
    axNonLin = subplot(2,1,2);
end
for timeIndex = 1:numel(timesteps)
    linState = [xLin; dxLin; theta1Lin; dtheta1Lin; theta2Lin; dtheta2Lin];
    nonLinState = [xNonLin; dxNonLin; theta1NonLin; dtheta1NonLin; theta2NonLin; dtheta2NonLin];
    FLin = K * linState;
    FNonLin = K * nonLinState;
    
    resultLin(timeIndex,:) = [timesteps(timeIndex), xLin, theta1Lin, theta2Lin, dxLin, dtheta1Lin, dtheta2Lin, FLin];
    resultNonLin(timeIndex,:) = [timesteps(timeIndex), xNonLin, theta1NonLin, theta2NonLin, dxNonLin, dtheta1NonLin, dtheta2NonLin, FNonLin];
    
    [xNonLin,dxNonLin,theta1NonLin,dtheta1NonLin,theta2NonLin,dtheta2NonLin] = simulateNonLinearSystem(xNonLin, dxNonLin, theta1NonLin, dtheta1NonLin, theta2NonLin, dtheta2NonLin, FNonLin, step, params);
    [xLin,dxLin,theta1Lin,dtheta1Lin,theta2Lin,dtheta2Lin] = simulateLinearSystem(xLin, dxLin, theta1Lin, dtheta1Lin, theta2Lin, dtheta2Lin, FLin, step, params);
    
    if GRAPHICAL_PLOT
        plotState(axLin, xLin, theta1Lin, theta2Lin, params.l1, params.l2);
        plotState(axNonLin, xNonLin, theta1NonLin, theta2NonLin, params.l1, params.l2);
    end
end
resultLin(end,:) = [timesteps(end) + step, xLin, theta1Lin, theta2Lin, dxLin, dtheta1Lin, dtheta2Lin, nan];
resultNonLin(end,:) = [timesteps(end) + step, xNonLin, theta1NonLin, theta2NonLin, dxNonLin, dtheta1NonLin, dtheta2NonLin, nan];
figure;
subplot 221;plot(resultLin(:,1),resultLin(:,2),'r',resultNonLin(:,1),resultNonLin(:,2),'b-.'); legend('Linear','Non-Linear'); title('x');
subplot 222;plot(resultLin(:,1),resultLin(:,3),'r',resultNonLin(:,1),resultNonLin(:,3),'b-.'); legend('Linear','Non-Linear'); title('t1');
subplot 223;plot(resultLin(:,1),resultLin(:,4),'r',resultNonLin(:,1),resultNonLin(:,4),'b-.'); legend('Linear','Non-Linear'); title('t2');
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