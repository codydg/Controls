%% Simuate Linear System
% Simulates the linear system for time "step" given power "u" with
% parameters "params" and initial values for x, theta1, theta2,
% and their respective time derivatives
function state_1 = simulateLinearSystem(state_0,F,step,params)
    % Input paramaters from state vector
    x_0 = state_0(1);
    dx_0 = state_0(2);
    theta1_0 = state_0(3);
    dtheta1_0 = state_0(4);
    theta2_0 = state_0(5);
    dtheta2_0 = state_0(6);
    
    % Input parameters from params struct
    g = params.g;
    m1 = params.m1;
    m2 = params.m2;
    M = params.M;
    l1 = params.l1;
    l2 = params.l2;
    
    % Calculate time-derivative of each input
    ddx = (-g*(m1*theta1_0 + m2*theta2_0) + F)/M;
    ddtheta1 = (-g*((M+m1)*theta1_0 + m2*theta2_0) + F)/(M*l1);
    ddtheta2 = (-g*(m1*theta1_0 + (M+m2)*theta2_0) + F)/(M*l2);
    
    % Increment inputs to calculate outputs
    dx_1 = dx_0 + ddx*step;
    dtheta1_1 = dtheta1_0 + ddtheta1*step;
    dtheta2_1 = dtheta2_0 + ddtheta2*step;
    x_1 = x_0 + dx_1*step;
    theta1_1 = theta1_0 + dtheta1_1*step;
    theta2_1 = theta2_0 + dtheta2_1*step;
    
    state_1 = [x_1; dx_1; theta1_1; dtheta1_1; theta2_1; dtheta2_1];
end