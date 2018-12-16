%% Run Simulation
%

%#ok<*UNRCH>
clc, clear, close all;

if ~exist('GRAPHICAL_PLOT','var')
    clc, clear, close all;
    
    GRAPHICAL_PLOT = false;
    PLOT_LINEAR = false;
    END_PLOT = true;
    MAKE_VIDEO = false;
    
    initState = [0;0;pi/16;0;0;0];

    Q = diag([1,1,10,1000,10,1000]);
    R = 0.00001;
    
    C = [1 0 0 0 0 0];
%     C = [1 0 0 0 1 0];
%     C = [1 0 1 0 1 0];
    
    finalX = 0;
%     finalX = 50;
end
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

stateLin = initState;
stateNonLin = initState;
estimatedStateLin = zeros(size(stateLin));
estimatedStateNonLin = estimatedStateLin;

AF = [0,1,0,0,0,0;0,0,-g*m1/M,0,-g*m2/M,0;0,0,0,1,0,0;0,0,-g*(M+m1)/(M*l1),0,-g*m2/(M*l1),0;0,0,0,0,0,1;0,0,-g*m1/(M*l2),0,-g*(M+m2)/(M*l2),0];
BF = [0;1/M;0;1/(M*l1);0;1/(M*l2)];

[K, S, e] = lqr(AF,BF,Q,R);

disp('Eigenvalues'' real parts from A_F - B_F * K');
fprintf('%f\n',real(e));

L = place(AF',C',real(e) * 10 + imag(e))';
eo = eig(AF - L*C);
disp(real(eo));

figure;
limits = [min(real(e)) - 0.1, 0.1, min(imag(e)) - 0.1,max(imag(e)) + 0.1];
plot(real(e),imag(e),'*',limits(1:2),[0,0],'k',[0,0],limits(3:4),'k');
xlabel('Real'); ylabel('Imaginary'); title('Poles of System');
axis(limits); grid on;
figure;
limits = [min(real(eo)) - 0.1, 0.1, min(imag(eo)) - 0.1,max(imag(eo)) + 0.1];
plot(real(eo),imag(eo),'*',limits(1:2),[0,0],'k',[0,0],limits(3:4),'k');
xlabel('Real'); ylabel('Imaginary'); title('Poles of System');
axis(limits); grid on;

%%
% Set up Time data
step = 0.0001; % Seconds
timesteps = 0:step:40-step;
resultLin = zeros(numel(timesteps) + 1, 14);
resultNonLin = resultLin;
if GRAPHICAL_PLOT
    fig = figure; fig.Position = [17 100 1200 728];
    ax1 = subplot(2,1,1); title(ax1,'Uncontrolled Lin');
    ax2 = subplot(2,1,2); title(ax2,'Uncontrolled NonLin');
end
if MAKE_VIDEO
    v = VideoWriter('video','MPEG-4');
    open(v);
end
FLin = zeros(numel(timesteps),1);
% FLin = 200 * ones(numel(timesteps),1);
FNonLin = FLin;
for timeIndex = 1:numel(timesteps)
    destimatedStateLin = AF*estimatedStateLin + BF*FLin(timeIndex) + L*C*(stateLin - estimatedStateLin);
    destimatedStateNonLin = AF*estimatedStateNonLin + BF*FNonLin(timeIndex) + L*C*(stateNonLin - estimatedStateNonLin);
    estimatedStateLin = estimatedStateLin + destimatedStateLin * step;
    estimatedStateNonLin = estimatedStateNonLin + destimatedStateNonLin * step;
    
    resultLin(timeIndex,:) = [timesteps(timeIndex), stateLin.', FLin(timeIndex), estimatedStateLin.'];
    resultNonLin(timeIndex,:) = [timesteps(timeIndex), stateNonLin.', FNonLin(timeIndex), estimatedStateNonLin.'];
    
    stateLin = simulateLinearSystem(stateLin, FLin(timeIndex), step, params);
    stateNonLin = simulateNonLinearSystem(stateNonLin, FNonLin(timeIndex), step, params);
    
    if GRAPHICAL_PLOT && mod(timeIndex, 10) == 0
        plotState(ax1, stateLin(1), stateLin(3), stateLin(5), params.l1, params.l2);
        plotState(ax2, stateNonLin(1), stateNonLin(3), stateNonLin(5), params.l1, params.l2);
        if MAKE_VIDEO
            writeVideo(v,getframe(gcf));
        end
    end
end
if MAKE_VIDEO
    close(v);
end
resultLin(end,:) = [timesteps(end) + step, stateLin.', nan, estimatedStateLin.'];
resultNonLin(end,:) = [timesteps(end) + step, stateNonLin.', nan, estimatedStateNonLin.'];
if (END_PLOT)
    figure('units','normalized','outerposition',[0 0 1 1]);
    N_SKIP = 10;
    subplot 321;plot(resultLin(1:end-N_SKIP,1),resultLin(1:end-N_SKIP,2) - resultLin(1:end-N_SKIP,9),'r',resultNonLin(1:end-N_SKIP,1),resultNonLin(1:end-N_SKIP,2) - resultNonLin(1:end-N_SKIP,9),'b-.'); legend('Linear','Non-Linear'); ylabel('X'); xlabel('Time'); title('Estimation Error X');
    subplot 323;plot(resultLin(1:end-N_SKIP,1),resultLin(1:end-N_SKIP,4) - resultLin(1:end-N_SKIP,11),'r',resultNonLin(1:end-N_SKIP,1),resultNonLin(1:end-N_SKIP,4) - resultNonLin(1:end-N_SKIP,11),'b-.'); legend('Linear','Non-Linear'); ylabel('Theta 1'); xlabel('Time'); title('Estimation Error Theta 1');
    subplot 325;plot(resultLin(1:end-N_SKIP,1),resultLin(1:end-N_SKIP,6) - resultLin(1:end-N_SKIP,13),'r',resultNonLin(1:end-N_SKIP,1),resultNonLin(1:end-N_SKIP,6) - resultNonLin(1:end-N_SKIP,13),'b-.'); legend('Linear','Non-Linear'); ylabel('Theta 2'); xlabel('Time'); title('Estimation Error Theta 2');
    subplot 322;plot(resultLin(1:end-N_SKIP,1),resultLin(1:end-N_SKIP,3) - resultLin(1:end-N_SKIP,10),'r',resultNonLin(1:end-N_SKIP,1),resultNonLin(1:end-N_SKIP,3) - resultNonLin(1:end-N_SKIP,10),'b-.'); legend('Linear','Non-Linear'); ylabel('d/dt X'); xlabel('Time'); title('Estimation Error d/dt X');
    subplot 324;plot(resultLin(1:end-N_SKIP,1),resultLin(1:end-N_SKIP,5) - resultLin(1:end-N_SKIP,12),'r',resultNonLin(1:end-N_SKIP,1),resultNonLin(1:end-N_SKIP,5) - resultNonLin(1:end-N_SKIP,12),'b-.'); legend('Linear','Non-Linear'); ylabel('d/dt Theta 1'); xlabel('Time'); title('Estimation Error d/dt Theta 1');
    subplot 326;plot(resultLin(1:end-N_SKIP,1),resultLin(1:end-N_SKIP,7) - resultLin(1:end-N_SKIP,14),'r',resultNonLin(1:end-N_SKIP,1),resultNonLin(1:end-N_SKIP,7) - resultNonLin(1:end-N_SKIP,14),'b-.'); legend('Linear','Non-Linear'); ylabel('d/dt Theta 2'); xlabel('Time'); title('Estimation Error d/dt Theta 2');
end
function plotState(ax, x, theta1, theta2, l1, l2)
    axes(ax);
    
    % Plot Ropes
    center1 = [x - l1*sin(theta1),-l1*cos(theta1)];
    center2 = [x - l2*sin(theta2),-l2*cos(theta2)];
    plot([center1(1), x, center2(1)], [center1(2), 0, center2(2)], 'k');
    hold on;
    
    % Plot Box
    boxWidth = 5; boxHeight = 5;
    rectangle('Position',[x-boxWidth/2, -boxHeight/2, boxWidth, boxHeight]);
    
    % Plot Circles
    massRadius = 2.5;
    rectangle('Position',[center1 - massRadius/2,massRadius,massRadius],'Curvature',1);
    rectangle('Position',[center2 - massRadius/2,massRadius,massRadius],'Curvature',1);
    
    % Adjust Axis
    xlim([-60+x,60+x]);
    ylim([-25,5]);
    drawnow; hold off;
end