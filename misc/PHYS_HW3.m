function [ output_args ] = PHYS_HW3( ~ )
%Matthew Howlett
%   Generate graph of the steady-state response, transient response, total
%   response and the force function with respect to time.

%declaration
t = 0:0.01:50;
m = 1.7;
k = 3.2;
b = 0.22;
F = 2.05;
omega_0 = 1.13;
A = -0.728;
B = -0.06247;

%calculation
gamma = b/m;
%omega_0 = sqrt(k/m);
omega_1 = sqrt(omega_0^2-(gamma/2)^2);
x_0 = (omega_1^2*(F/m))/(sqrt((omega_0^2-omega_1^2)^2+gamma^2*omega_1^2));
phi = atan((omega_1*gamma)/(omega_0^2-gamma^2));
x_1 = x_0*cos(omega_1*t-phi); 
x_2 = (A*exp((-gamma/2)*t).*cos(omega_1*t))+(B*exp((-gamma/2)*t).*sin(omega_1*t));
x_t = x_1+x_2;
f_t = F*cos(omega_1*t);
p = (F^2)/(2*m*gamma*((4*(omega_1-omega_0)^2)/(gamma^2+1)));



%print to screen
disp('Damping coefficient:');
disp(gamma);
disp('frequency:');
disp(omega_1);
disp('Steady State initial displacement:');
disp(x_0);
disp('Phase angle:');
disp(phi);
disp('Power dissipated:')
disp(p);
disp('Unkown A:');
disp(A);
disp('Unknown B:')
disp(B);

%plot
plot(t, x_1, '-r','LineWidth',1); hold on;
plot(t, x_2, '-b','LineWidth',1);
plot(t, x_t, '-g','LineWidth',1);
plot(t, f_t, '-m','LineWidth',1);
hold off;


end

