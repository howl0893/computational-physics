%% Euler_Cromer_OscillatoryMotion.m
%
% algorithm designed to plot position and velocity for a harmonic
% oscillator problem using the Last Point Approximation method (LPA).

%% Define Parameters
iterations = 1000;      % number of iterations in LPA 
m = 0.1;                % mass of object in motion [kg]
x0 = 0.1;               % initial position is at equilibirum
v0 = 0;                 % mass begins motion at rest
dt = 0.001;             % step time inbetween point approximations
k = 10;                 % Spring constant [N/m]
T = (2*pi)*(m/k).^0.5;  % theoretical period
t = 0:((5*T)/1000):5*T; % time elapsed
x = ones(1,iterations); % intiialize position vector to store values
v = ones(1,iterations); % initialize velocity vector to store values

%% Calculations

x(1) = x0*x(1);             % set initial position
v(1) = v0*v(1);             % set initial velocity

tic                         % start stopwatch to measure computation time

for i = (1:1:iterations)
    v(i+1) = v(i) - (k/m)*x(i).*dt; % Approximate velocity
    x(i+1) = x(i)+v(i+1).*dt;       % Approximate position 
end

toc

%% Plotting
figure(1);                             % create first figure
hold on;

title('Position of Mass versus Time'); % Label title
ylabel('Position [m]');                % Label y axis
xlabel('Time [s]');                    % Label x axis
xlim([0 (5*T)]);                       % set x limits

plot(t,x);                             %plot position vs time
plot(t, time);

hold off;

figure(2);                             % create second figure
hold on;

title('Velocity of Mass versus Time'); % Label title
ylabel('Velocity [m/s]');              % Label y axis
xlabel('Time [s]');                    % Label x axis
xlim([0 (5*T)]);                       % set x limits

plot(t,v);                             %plot velocity vs time
plot(time);

hold off;
