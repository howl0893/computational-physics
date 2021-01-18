%% BallisticTrajectory
%
% Comparing the trajectory for a ballistic projectile 
% using kinematics and Euler's forward method

%% Declaration of Parameters
dt = 0.1;            % step size
Ax = 0;              % x Acceleration
Ay= 9.8;             % y Acceleration
Vx = 50*ones(1,110); % x initial velocity
Vy = 50*ones(1,110); % y initial velocity
t = (1:1:111);       % time of flight
X = ones(1,110);     % X vector for theoretical solution
Y = ones(1,110);     % Y vector for theoretical solution

%For kinematics calculation
t1 = (1:1:11);       % time of flight
Vx1 = 50*ones(1,11); % x initial velocity
Vy1 = 50*ones(1,11); % y initial velocity

%% Kinematics Calculated Solution
X1 = (Vx1.*t1)+((1/2)*(Ax)*t1.^2);   % displacement calculation
Y1 = (Vy1.*t1)+((1/2)*(-Ay).*t1.^2); % height calculation
X1(1) = 0;                           % set initial X to be zero
Y1(1) = 0;                           % set initial Y to be zero


%% Euler Forward Method Solution
for i = (1:1:110)
    X(i+1) = X(i)+ Vx(i).*dt;   % displacement
    Vx(i+1) = Vx(i);            % velocity in X direction
    Y(i+1) = Y(i)+ Vy(i).*dt;   % height
    Vy(i+1) = Vy(i) - Ay.*dt;   % velocity in Y direction
end

%% Plotting

% make first figure
figure(1);                                  
hold on;

%Plot finite differences vs time
plot(t,X,'r');
plot(t,Vx,'b');
plot(t,Y,'g');
plot(t,Vy,'y');

%Label axes and title
xlabel('Time [s]');
ylabel('Height [m]');
title('Finite Difference Equation Against Time');

%Inlcude legend
legend('position X', 'velocity X', 'position Y', 'velocity Y');
hold off;

%Make second figure
figure(2); 
hold on;

plot(X,Y, 'r-'); %Plot theoretical

plot(X1,Y1, 'b-'); 

%Label axes and title
xlabel('Distance [m]'); 
ylabel('Height [m]');
title('Theoretical vs. Calculated Trajectory for Ballistic Projectile');

%Include legend
legend('Theoretical', 'Calculated');

hold off;

