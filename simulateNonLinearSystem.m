%% Simuate NonLinear System
% Simulates the nonlinear system for time "step" given power "u" with
% parameters "params" and initial values for x, theta1, theta2,
% and their respective time derivatives
function state_1 = simulateNonLinearSystem(state_0,F,step,params)
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
    
    J1 = m1*sin(theta1_0)*(l1*dtheta1_0^2 + g*cos(theta1_0));
    J2 = m2*sin(theta2_0)*(l2*dtheta2_0^2 + g*cos(theta2_0));
    
    % Calculate time-derivative of each input
    ddx = (F - J1 - J2)/(M + m1*sin(theta1_0)^2 + m2*sin(theta2_0)^2);
    ddtheta1 = -(g*sin(theta1_0)/l1) + (cos(theta1_0)/l1)*ddx;
    ddtheta2 = -(g*sin(theta2_0)/l2) + (cos(theta2_0)/l2)*ddx;
    
    % Increment inputs to calculate outputs
    dx_1 = dx_0 + ddx*step;
    dtheta1_1 = dtheta1_0 + ddtheta1*step;
    dtheta2_1 = dtheta2_0 + ddtheta2*step;
    x_1 = x_0 + dx_1*step;
    theta1_1 = theta1_0 + dtheta1_1*step;
    theta2_1 = theta2_0 + dtheta2_1*step;
    
    state_1 = [x_1; dx_1; theta1_1; dtheta1_1; theta2_1; dtheta2_1];
end