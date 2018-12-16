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
    
    initState = [0;10;0;pi/8;0;-pi/8];    
%     initState = [0;0;0;0;0;0];

    Q = diag([1,1,10,1000,10,1000]);
    R = 0.00001;
    
%     C = [1 0 0 0 0 0];
%     C = [1 0 0 0 1 0];
    C = [1 0 1 0 1 0];
    
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
controlledStateLin = initState;
controlledStateNonLin = initState;
estimatedStateLin = diag(C)*initState;
estimatedStateNonLin = diag(C)*initState;
estimatedConStateLin = diag(C)*initState;
estimatedConStateNonLin = diag(C)*initState;

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
step = 0.01; % Seconds
timesteps = 0:step:40-step;
resultLin = zeros(numel(timesteps) + 1, 14);
resultNonLin = resultLin;
controlledResultLin = resultLin;
controlledResultNonLin = resultLin;
if GRAPHICAL_PLOT
    fig = figure; fig.Position = [17 100 1200 728];
    if PLOT_LINEAR
        ax1 = subplot(2,1,1); title(ax1,'Uncontrolled Lin');
        ax2 = subplot(2,1,2); title(ax2,'Controlled Lin');
    else
        ax1 = subplot(2,1,1); title(ax1,'Uncontrolled NonLin');
        ax2 = subplot(2,1,2); title(ax2,'Controlled NonLin');
    end
end
if MAKE_VIDEO
    v = VideoWriter('video','MPEG-4');
    open(v);
end
FLin = 0;
FNonLin = 0;
FConLin = 0;
FConNonLin = 0;
for timeIndex = 1:numel(timesteps)
    destimatedStateLin = AF*estimatedStateLin + BF*FLin + diag(L)*diag(C)*(stateLin - estimatedStateLin);
    destimatedStateNonLin = AF*estimatedStateNonLin + BF*FNonLin + diag(L)*diag(C)*(stateNonLin - estimatedStateNonLin);
    destimatedConStateLin = AF*estimatedConStateLin + BF*FConLin + diag(L)*diag(C)*(controlledStateLin - estimatedConStateLin);
    destimatedConStateNonLin = AF*estimatedConStateNonLin + BF*FConNonLin + diag(L)*diag(C)*(controlledStateNonLin - estimatedConStateNonLin);
    estimatedStateLin = estimatedStateLin + destimatedStateLin * step;
    estimatedStateNonLin = estimatedStateNonLin + destimatedStateNonLin * step;
    estimatedConStateLin = estimatedConStateLin + destimatedConStateLin * step;
    estimatedConStateNonLin = estimatedConStateNonLin + destimatedConStateNonLin * step;
    FLin = 0;
    FNonLin = 0;
    FConLin = -K * (estimatedConStateLin - [finalX;zeros(5,1)]);
    FConNonLin = -K * (estimatedConStateNonLin - [finalX;zeros(5,1)]);
    
    resultLin(timeIndex,:) = [timesteps(timeIndex), stateLin.', FLin, estimatedStateLin.'];
    resultNonLin(timeIndex,:) = [timesteps(timeIndex), stateNonLin.', FNonLin, estimatedStateNonLin.'];
    controlledResultLin(timeIndex,:) = [timesteps(timeIndex), controlledStateLin.', FConLin, estimatedConStateLin.'];
    controlledResultNonLin(timeIndex,:) = [timesteps(timeIndex), controlledStateNonLin.', FConNonLin, estimatedConStateNonLin.'];
    
    stateLin = simulateLinearSystem(stateLin, FLin, step, params);
    stateNonLin = simulateNonLinearSystem(stateNonLin, FNonLin, step, params);
    controlledStateLin = simulateLinearSystem(controlledStateLin, FConLin, step, params);
    controlledStateNonLin = simulateNonLinearSystem(controlledStateNonLin, FConNonLin, step, params);
    
    if GRAPHICAL_PLOT && mod(timeIndex, 10) == 0
        if PLOT_LINEAR
            plotState(ax1, stateLin(1), stateLin(3), stateLin(5), params.l1, params.l2);
            plotState(ax2, controlledStateLin(1), controlledStateLin(3), controlledStateLin(5), params.l1, params.l2);
        else
            plotState(ax1, stateNonLin(1), stateNonLin(3), stateNonLin(5), params.l1, params.l2);
            plotState(ax2, controlledStateNonLin(1), controlledStateNonLin(3), controlledStateNonLin(5), params.l1, params.l2);
        end
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
controlledResultLin(end,:) = [timesteps(end) + step, controlledStateLin.', nan, estimatedConStateLin.'];
controlledResultNonLin(end,:) = [timesteps(end) + step, controlledStateNonLin.', nan, estimatedConStateNonLin.'];
if (END_PLOT)
    figure('units','normalized','outerposition',[0 0 1 1]);
    subplot 321;plot(resultLin(:,1),resultLin(:,2),'r',resultNonLin(:,1),resultNonLin(:,2),'b-.'); legend('Linear','Non-Linear'); ylabel('X'); xlabel('Time'); title('Uncontrolled X');
    subplot 323;plot(resultLin(:,1),resultLin(:,4),'r',resultNonLin(:,1),resultNonLin(:,4),'b-.'); legend('Linear','Non-Linear'); ylabel('Theta 1'); xlabel('Time'); title('Uncontrolled Theta 1');
    subplot 325;plot(resultLin(:,1),resultLin(:,6),'r',resultNonLin(:,1),resultNonLin(:,6),'b-.'); legend('Linear','Non-Linear'); ylabel('Theta 2'); xlabel('Time'); title('Uncontrolled Theta 2');
    subplot 322;plot(controlledResultLin(:,1),controlledResultLin(:,2),'r',controlledResultNonLin(:,1),controlledResultNonLin(:,2),'b-.'); legend('Linear','Non-Linear'); ylabel('X'); xlabel('Time'); title('Controlled X');
    subplot 324;plot(controlledResultLin(:,1),controlledResultLin(:,4),'r',controlledResultNonLin(:,1),controlledResultNonLin(:,4),'b-.'); legend('Linear','Non-Linear'); ylabel('Theta 1'); xlabel('Time'); title('Controlled Theta 1');
    subplot 326;plot(controlledResultLin(:,1),controlledResultLin(:,6),'r',controlledResultNonLin(:,1),controlledResultNonLin(:,6),'b-.'); legend('Linear','Non-Linear'); ylabel('Theta 2'); xlabel('Time'); title('Controlled Theta 2');
    figure; plot(controlledResultLin(:,1),controlledResultLin(:,8),'r',controlledResultNonLin(:,1),controlledResultNonLin(:,8),'b-.'); legend('Linear','Non-Linear'); ylabel('Control Input'); xlabel('Time'); title('Control Input (Force)');
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