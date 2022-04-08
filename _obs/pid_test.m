% Create msd transfer function
m = .001; % kg
k = 1; % N/m
c = 0.0001; % N*s/m

s = tf('s');
P = 1/(m*s^2 + k*s + c);

step(P)

figure
Kp = 300;
Ki = 100;
Kd = 2;
C = pid(Kp,Ki,Kd)
T = feedback(C*P,1);

t = 0:0.0001:0.05;
bode(T)
%step(T,t)
hold off
%pidTuner(P,C)

