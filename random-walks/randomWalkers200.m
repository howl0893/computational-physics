%% Random Walker
%
% This algorithm is designed to plot the random walk of six "walker". 

%% Parameters
x(1) = 0;         % initial x position
y(1) = 0;         % initial y position
maxStep = 500;    % max step size [nm]
maxIndex = 50;    % max iterations
timeStep = 0.001; % time step in seconds
walkerNum = 200;  % Number of walkers
xMaster = [];     % used to store all the walkers x-coordinates
yMaster = [];     % used to store all the walkers y-coordinates
time = 0:timeStep:(0.001*50);

%% Method
for walker = 1:1:walkerNum
    for i = 1:maxIndex
        first = rand;                           % first random number
        second = rand;                          % second random number
        stepsize = first*maxStep;               % step size of random walk
        angle = second * 2*pi;                  % direction of random walk
        x(i+1) = x(i) + stepsize * cos(angle);  % next x position
        y(i+1) = y(i) + stepsize * sin(angle);  % next y position
    end
    xMaster(walker,:) = [walker, x];    % store all walkers x-coordinate
    yMaster(walker,:) = [walker, y];    % store all walkers y-coordinate
end

for tIndex = 2:1:length(time)
    for i = 2:walkerNum
        rSquared(tIndex) = xMaster(i,tIndex).^2+yMaster(i,tIndex).^2;
    end
    rSquaredAvg(tIndex) = sum(rSquared)/walkerNum;
end

%% Plotting
figure(1);
hold on;
for i = 1:walkerNum
    plot(xMaster(i,2:52),yMaster(i,2:52),'-o');      % plot
end
xlabel('X [nm]');        % label
ylabel('Y [nm]');        % label
title('Random Walker');  % label

figure(2);
hold on;
plot(time,rSquaredAvg);
title('MSD of all the Walkers');
xlabel('Time [s]');
ylabel('Displacement [nm]');