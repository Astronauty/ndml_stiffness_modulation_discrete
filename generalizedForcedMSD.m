%% Model of Generalized 1D Discrete Nonreciprocal Vibration Transmission
% Daniel Nguyen
% 4/11/20


% Global variable declarations 
% TODO Change to passed parameters once methods have been finalized
%clear all
close all
%clc

global N d 
global k_wavenumber k_angularfreq c_wavenumber c_angularfreq A_k A_c

N = 600; % Number of discrete masses/coordinates in the system
d = .1; % Distance between masses

tic
%% Set stiffness and damping modulation rates
A_k = 5;
k_wavenumber = 0.1*(2*pi);
k_angularfreq = 0.1*(2*pi);

A_c = 0;
c_wavenumber = 0*(2*pi);
c_angularfreq = 0*(2*pi);

%% State Variable Integration
y0 = [0; zeros(N-1,1); 1; zeros(N-1,1)]; % Initial conditions, first half of vector is initial position and second half is initial velocity

ts = [0 20]; % Time domain (s)

%[t,y] = ode45('f',ts,y0);
 options = odeset('AbsTol',1e-3,'RelTol',1e-3,'Stats','on');
[t,y] = ode45('f',ts,y0,options);


%% Plotting
tiledlayout(6,2); % Requires R2019b or later

% Plot Initial k and c
x = [1:1:N]*d;
[c,C] = get_damping(0,k_wavenumber,k_angularfreq);
[k,K] = get_stiffness(0,k_wavenumber,k_angularfreq);

nexttile
scatter(x,k);
title('Initial Stiffness Space-Dependence')
xlabel('Distance (m)') 
ylabel('Stiffness (N/m)') 

nexttile
scatter(x,c);
title('Initial Damping Space-Dependence')
xlabel('Distance (m)') 
ylabel('Damping (kg/s)') 

% Plot Displacement
nexttile([1 2])
plot(t,y(:,1),t,y(:,25));
axis on
title('Displacement vs. Time of Discrete Masses')
xlabel('Time (s)') 
ylabel('Displacement (m)') 
legend({'Mass 1','Mass 25'},'Location','northeast')


% Plot Velocity
nexttile([1 2])
plot(t,y(:,N+1),t,y(:,2*N));
title('Velocity vs. Time of Discrete Masses')
xlabel('Time (s)') 
ylabel('Velocity (m/s)') 
legend({'Mass 1','Mass 25'},'Location','northeast')
axis on

%% FFT Analysis
nexttile([1 2])

% Store displacement data
m1_disp = y(:,1);
m2_disp = y(:,25);


% Create filtering window between x = filter_start and x = filter_end
filter_start = 80;
filter_end = 260;

w = [zeros(filter_start,1); window(@tukeywin,filter_end-filter_start); zeros(length(t)-filter_end,1)];
m2_disp_filtered = y(:,25).*w;

plot(t,w,'r');
hold on
plot(t,m2_disp,'c');
hold on
plot(t,m2_disp_filtered,'b');

title('Windowing Function')
xlabel('Time (s)') 
ylabel('Displacement (m)') 

legend({'Windowing Function','Mass 25 Disp','Mass 25 Disp, Windowed'},'Location','northeast')


nexttile([1 2])
[f1,Sig1] = FFT_perso(t,m1_disp); % Perform FFT on 1st mass
[f2,Sig2] = FFT_perso(t,m2_disp); % Perform FFT on 25th mass
[f2_filt,Sig2_filt] = FFT_perso(t,m2_disp_filtered); % Perform FFT on 25th mass



plot(f1,abs(Sig1),'r')
hold on
plot(f2,abs(Sig2),'c')
hold on
plot(f2_filt,abs(Sig2_filt),'b')

legend({'Mass 1 (Raw)','Mass 25 (Raw)','Mass 25 (Windowed)'},'Location','northeast')
xlim([0,4]);
title('Displacement FFT')
xlabel('Frequency (Hz)') 
ylabel('Amplitude') 



%% Plot total system energy
    % M = eye(N);
    % E_displacement = 0.5*K.*(y(500,1:N).^2);
    % E_velocity = 0.5*M.*(y(500,N+1:2*N));

[E_displacement, E_velocity, E_total] = getTotalEnergy(t,y,k);
% 
% nexttile
% plot(t,E_displacement);
% title('System Potential Energy')
% xlabel('Time (s)') 
% ylabel('U (J)') 
% 
% nexttile
% plot(t,E_velocity)
% title('System Kinetic Energy')
% xlabel('Time (s)') 
% ylabel('T (J)') 

nexttile([1 2])
plot(t,E_total)
title('System Total Energy')
xlabel('Time (s)') 
ylabel('Total Energy(J)') 

toc

%% Forward vs Reverse
figure
loglog(f_forward,sig_forward,'b')
hold on
loglog(f_reverse,sig_reverse,'r')
title('Forward/Reverse FFT')
xlabel('Frequency (Hz)') 
ylabel('Amplitude') 

legend({'Mass 25 Forward','Mass 25 Reverse'},'Location','northeast')

