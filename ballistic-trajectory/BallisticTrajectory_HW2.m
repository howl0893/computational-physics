%% BallisticTrajectory_HW2
%
% This algorithm is designed to plot the trajectory and velocity
% for a cannoball projectile of both the kinematic calculation
% and the theoretical calculation based on Euler's Forward method
% The commented out portion of the code is designed to plot the 
% trajectory of a cannonball at various angles.

%% Declaration of Parameters

% variables for kinematics calculation
x0=0;               % initial x posiition
y0=0;               % initial y position
angle = (pi./4);    % angle of 45 degrees
angle2 = (pi./6);   % angle of 30 degrees
angle3 = (pi./3);   % angle of 60 degrees
v = 50;             % initial velocity
g = 9.81;           % acceleration due to gravity
ax = 0;             % acceleration in the x direction
t = 0:0.1:20;       % time elapsed

% variables for theoretical calculation
dt = .009;                       % step size
Ax = 0;                           % x Acceleration
Ay= 9.8;                          % y Acceleration
Vx = 50*ones(1,1000)*cos(angle3); % x initial velocity
Vy = 50*ones(1,1000)*sin(angle3); % y initial velocity
t1 = (1:1:1001);                  % time of flight
X = ones(1,1000);                 % X vector for theoretical solution
Y = ones(1,1000);                 % Y vector for theoretical solution

%% Kinematic Calculations
x = x0+v*cos(angle)*t+(ax*t.^2)/2;   % x calculation at 45 degrees
y = y0+v*sin(angle)*t-(g*t.^2)/2;    % y calculation at 45 degrees
x2 = x0+v*cos(angle2)*t+(ax*t.^2)/2; % x calculation at 30 degrees
y2 = y0+v*sin(angle2)*t-(g*t.^2)/2;  % y calculation at 30 degrees
x3 = x0+v*cos(angle3)*t+(ax*t.^2)/2; % x calculation at 60 degrees
y3 = y0+v*sin(angle3)*t-(g*t.^2)/2;  % y calculation at 60 degrees

%% Euler Forward Method Solution
for i = (1:1:1000)
    X(i+1) = X(i)+ Vx(i).*dt;   % displacement
    Vx(i+1) = Vx(i);            % velocity in X direction
    Y(i+1) = Y(i)+ Vy(i).*dt;   % height
    Vy(i+1) = Vy(i) - Ay.*dt;   % velocity in Y direction
end

%% Plotting

%plot trajectory 
figure(1)
hold on;

plot(x,y);          % plot trajectory at 45 degrees
plot(x2,y2);        % plot trajectory at 30 degrees
plot(x3,y3);        % plot trajectory at 60 degrees

plot(X,Y);         % plot theoretical

ylim([0 100]);          % set y limits
xlim([0 250]);          % set x limits

xlabel('Distance [m]'); % label x axis
ylabel('Height [m]');   % label y axis

title('Cannonball Trajectory - Theoretical vs. Calculated');       % label title
%legend('Kinematic', 'Theroretical');                               % create legend
legend('45 degrees', '30 degrees', '60 degrees', 'theroretical'); % create legend

hold off;

% plot velocity
% figure(2)
% hold on;
% 
% plot(t1, Vx,'k'); % plot x velocity
% plot(t1, Vy,'m'); % plot y velocity
% 
% ylim([0 100]);    % set y limits
% xlim([0 500]);    % set x limits
% 
% xlabel('Time [s]'); % label x axis
% ylabel('Velocity [m/s]');   % label y axis
% 
% title('Cannonball Velocity vs. Time');       % label title
% 
% legend('X-Velocity', 'Y-Velocity'); % Create Legend
% 
% hold off;

