function [ output_args ] = PHYS_HW5( ~ )
%Matthew Howlett

%declaration
x = -10:0.1:10;
t = 0;
T = 100;
v = 25;
b = 0.05;


%calculation
y = b.^3/b.^2+(x-v*t).^2;


%print to screen
%disp('Transverse Wave Traveling along a String:');
%disp(y);


%plot
plot(x, y, '-r','LineWidth',1);


end
