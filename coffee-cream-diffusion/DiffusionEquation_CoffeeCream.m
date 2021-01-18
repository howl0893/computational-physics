%% DiffusionEquation_CoffeeCream.m
%
% This algorithm is designed to create a simulation of a drop of cream
% diffusing into a cup of coffee over time - limited to two dimensions.

%% Parameters
L=15;                            % length of side of mug, cm
nmax = 1000;                     % time step maximum
Pinitial=100;                    % initial concentration of creamer, g/cm^2
Pold = zeros(L,L,1);             % old concentration
Pchange=zeros(L,L,1);            % change in concentration
Pfigure = zeros(L,L,1);          % concentration to be plotted
D = 0.25;                        % diffusion constant
deltat = 0.25;                   % change in time
deltax = 1;                      % change in x
Pold((L+1)/2,(L+1)/2,1)=Pinitial;% put the initial concentration in the middle of the mug


%% Calculations
% Note: Press the space bar to step through the diffusion process - hold
% down the space bar to speed up the process. ctrl+c to exit.
for n=1:nmax
    for i = 2:1:14
        for j = 2:1:14
            % finite difference equation
            Pchange(i,j,n+1) = Pold(i,j,n)+(deltat*D/deltax.^2).*(Pold(i+1,j,n)+Pold(i-1,j,n)+Pold(i,j+1,n)+Pold(i,j-1,n)-4.*Pold(i,j,n));
        end
    end
    
    %Plotting
    Pold = Pchange;          % set Pold
    figure(1);               % create figure 
    Pfigure = Pold(:,:,n+1); 
    surf(Pfigure);           % plot figure
    axis([0 15 0 15 0 100]); 
    colorbar;
    grid on;
    title('2D Gaussian Curve of expected Cream Diffusion'); % title plot
    xlabel('X direction (cm)');                             % label x axis
    ylabel('Y direction (cm)');                             % label y axis
    zlabel('Concentration (g/cm^2)');                       % label z axis% set axis
    disp('Press space to advance')
    pause                    % pause each iteration
end 
    
sum(sum(Pold));