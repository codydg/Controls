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
    
%     initState = [0;1;0;pi/16;0;-pi/16];
    initState = [0;0;0;0;0;0];
    
%     L_Goal = [-20;-15;-10;-5;-3.5;-2]*0.1;
    L_Goal = [-20;-19;-10;-9;-5;-4]*0.15;
%     L_Goal = linspace(-0.3,-0.15,6)*1;

    step = 0.0001; % Seconds
    timesteps = 0:step:50-step;
    F_Summary = zeros(numel(timesteps),1);
    F_Summary(20/step:end) = 100;
    lim = 0.04*step*F_Summary(end);
    
    C = [1 0 0 0 0 0];
%     C = [1 0 0 0 1 0];
%     C = [1 0 1 0 1 0];
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
estimatedStateLin = diag(C) * initState;
estimatedStateNonLin = estimatedStateLin;

AF = [0,1,0,0,0,0;0,0,-g*m1/M,0,-g*m2/M,0;0,0,0,1,0,0;0,0,-g*(M+m1)/(M*l1),0,-g*m2/(M*l1),0;0,0,0,0,0,1;0,0,-g*m1/(M*l2),0,-g*(M+m2)/(M*l2),0];
BF = [0;1/M;0;1/(M*l1);0;1/(M*l2)];

L = place(AF',C',L_Goal).';
eo = eig(AF - L*C);
disp(real(eo));

figure;
limits = [min(real(eo)) - 0.1, 0.1, min(imag(eo)) - 0.1,max(imag(eo)) + 0.1];
plot(real(eo),imag(eo),'*',limits(1:2),[0,0],'k',[0,0],limits(3:4),'k');
xlabel('Real'); ylabel('Imaginary'); title('Poles of System');
axis(limits); grid on;

% Set up Time data
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
FLin = F_Summary;
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
    
    if GRAPHICAL_PLOT && mod(timeIndex, 50) == 0
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
    subplot 321;plot(resultLin(1:end,1),resultLin(1:end,2) - resultLin(1:end,9),'r',resultNonLin(1:end,1),resultNonLin(1:end,2) - resultNonLin(1:end,9),'b-.'); legend('Linear','Non-Linear'); ylabel('X'); xlabel('Time'); title('Estimation Error X'); 
    if ~isnan(lim)
        ylim([-lim,lim]);
    end
    subplot 323;plot(resultLin(1:end,1),resultLin(1:end,4) - resultLin(1:end,11),'r',resultNonLin(1:end,1),resultNonLin(1:end,4) - resultNonLin(1:end,11),'b-.'); legend('Linear','Non-Linear'); ylabel('Theta 1'); xlabel('Time'); title('Estimation Error Theta 1');
    if ~isnan(lim)
        ylim([-lim,lim]/20);
    end
    subplot 325;plot(resultLin(1:end,1),resultLin(1:end,6) - resultLin(1:end,13),'r',resultNonLin(1:end,1),resultNonLin(1:end,6) - resultNonLin(1:end,13),'b-.'); legend('Linear','Non-Linear'); ylabel('Theta 2'); xlabel('Time'); title('Estimation Error Theta 2');
    if ~isnan(lim)
        ylim([-lim,lim]/20);
    end
    subplot 322;plot(resultLin(1:end,1),resultLin(1:end,3) - resultLin(1:end,10),'r',resultNonLin(1:end,1),resultNonLin(1:end,3) - resultNonLin(1:end,10),'b-.'); legend('Linear','Non-Linear'); ylabel('d/dt X'); xlabel('Time'); title('Estimation Error d/dt X');
    if ~isnan(lim)
        ylim([-lim,lim]/10);
    end
    subplot 324;plot(resultLin(1:end,1),resultLin(1:end,5) - resultLin(1:end,12),'r',resultNonLin(1:end,1),resultNonLin(1:end,5) - resultNonLin(1:end,12),'b-.'); legend('Linear','Non-Linear'); ylabel('d/dt Theta 1'); xlabel('Time'); title('Estimation Error d/dt Theta 1');
    if ~isnan(lim)
        ylim([-lim,lim]/20);
    end
    subplot 326;plot(resultLin(1:end,1),resultLin(1:end,7) - resultLin(1:end,14),'r',resultNonLin(1:end,1),resultNonLin(1:end,7) - resultNonLin(1:end,14),'b-.'); legend('Linear','Non-Linear'); ylabel('d/dt Theta 2'); xlabel('Time'); title('Estimation Error d/dt Theta 2');
    if ~isnan(lim)
        ylim([-lim,lim]/20);
    end
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