%% Unit 5 Project - Self-Avoiding Random Walk
%
% This program is designed to model a self-avoiding random walk on a
% two-dimensional lattice.

%% Parameters
numWalkers = 200; % number of walkers
dx = 1;           % change in x
dy = 1;           % change in y
stepSize = 1;     % step size
maxSteps = 50;   % maximum number of steps
dt = 0.01;        % change in time
time = 0;         % time

xMax = maxSteps;  % x maximum
yMax = maxSteps;  % y maximum

x = zeros(numWalkers,xMax); % initialize x-matrix
y = zeros(numWalkers,yMax); % initialize y-matrix
stuckWalkers = 0;           % number of 'stuck' walkers
i = 1;                      % used for indexing

%% Method
for row = 1:numWalkers
    
    % space for each walker
    stopWalker = false;                       % flag to stop walker - initially set to false
    space = zeros(maxSteps*2-1,maxSteps*2-1); % initialize total area walkers can walk
    space(maxSteps,maxSteps) = 1;             % set the center of space matrix to 1
    
    for column = 1:stepSize:maxSteps
        
       stopWalker = false;  % flag to stop walker - initially set to false
       hasMoved = 0;        % flag to indicate walker has moved
       iteration = 0;       % initial iteration number
       
       % looks in each direction surrounding one point
       % 1 if another walker is there, 0 if not
       left = space(maxSteps+y(row,column), maxSteps+x(row,column)-stepSize);
       right = space(maxSteps+y(row,column), maxSteps+x(row,column)+stepSize);
       up = space(maxSteps+y(row,column)-stepSize, maxSteps+x(row,column));
       down = space(maxSteps+y(row,column)+stepSize, maxSteps+x(row,column));
       
       while hasMoved == 0
           
           random = ceil(4*rand()); % create random number - round up to whole integer
           
           % random number is mulitiplied by a factor of 4 for each
           % direction it can move. Based on this random value, the particle will move... 
           
           % move left
           if random == 1 && left == 0
               y(row,column+1) = y(row,column);
               x(row,column+1) = x(row,column) - stepSize;
               space(maxSteps+y(row,column), maxSteps+x(row,column)-1) = 1;
               hasMoved=1;
           % move right
           elseif random == 2 && right == 0
               y(row,column+1) = y(row,column);
               x(row,column+1) = x(row,column) + stepSize;
               space(maxSteps+y(row,column), maxSteps+x(row,column)+1) = 1;
               hasMoved=1;
           % move up
           elseif random == 3 && up == 0
               y(row,column+1) = y(row,column) - stepSize;
               x(row,column+1) = x(row,column);
               space(maxSteps+y(row,column)-1, maxSteps+x(row,column)) = 1;
               hasMoved = 1;
           % move down    
           elseif random == 4 && down == 0
               y(row,column+1) = y(row,column) + stepSize;
               x(row,column+1) = x(row,column);
               space(maxSteps+y(row,column)+1, maxSteps+x(row,column)) = 1;
               hasMoved = 1;
           end
           
           % break while loop after 50 iterations
           iteration = iteration + 1;
           if iteration >= 50
               stopWalker = true;
               break
           end
           
       end
       
       rSquared(row,column+1) = (x(row,column+1)-x(row,1))^2 + (y(row,column+1)-y(row,1))^2; %.
       time(column+1) = time(column) + dt;
       
       % rebreak
       if stopWalker == true
           stuckWalkers(1,i) = row;
           i = i + 1;
           break
       end
       
    end
    
    % plot walker motion 
    if stopWalker == true
        figure(1);
        plot(x(row,1:column),y(row,1:column), '-');
        hold on;
    else
        figure(1);
        plot(x(row,1:maxSteps),y(row,1:maxSteps), '-');
        hold on;
    end
end

% calculate <r^2>, D, and alpha values
rSquaredAvg = sum(rSquared)./numWalkers;
D = rSquaredAvg./(4*time); 
alpha = (log(max(rSquaredAvg)./max(D)))./(log(time(maxSteps+1))); 
% rSquaredFit = D.*time.^alpha;

% display D and alpha
disp('Diffusion Coefficent = ');
disp(D);
disp('Alpha = ');
disp(alpha);

%% Plotting
% figure 1 settings
grid on;
grid minor;
titleStr = [num2str(numWalkers), ' Self-Avoiding Random Walker(s)'];
title(titleStr);
xlabel('X-Position [m]');
ylabel('Y-Position [m]');

figure(2);
plot(time,rSquaredAvg);
% plot(time,rSquaredFit);
xlim([0 time(maxSteps+1)])
titleStr2 = ['Mean Squared Displacement vs. Time of ', num2str(numWalkers), ' Walker(s)'];
title(titleStr2);
xlabel('Time [s]');
ylabel('Mean Squared Displacement [m^2]');

%% Fit: 'Mean Squared Displacement'.
[xData, yData] = prepareCurveData( time, rSquaredAvg );

% Set up fittype and options.
ft = fittype( 'poly1' );

% Fit model to data.
[fitresult, gof] = fit( xData, yData, ft );

% Plot fit with data.
figure(3);
h = plot( fitresult, xData, yData );
legend( h, '<r^2>', 'Linear Fit', 'Location', 'NorthEast' );

%set axe limit
xlim([0 0.5]);

% Label axes
title('Mean Squared Displacement vs. time with curve fit');
xlabel('time [s]');
ylabel('Displacement [m^2]');
grid on


