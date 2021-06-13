%% Model of Generalized 1D Discrete Nonreciprocal Vibration Transmission
% Daniel Nguyen
% 4/11/20
%% Initialize default vars
tic
N = 400; % Number of masses
d = 0.1; % Distance between masses
y0 = [0; zeros(N-1,1); 0; zeros(N-1,1)]; % Initial conditions (currently 0 position, 0 velocity)
ts = [0 20]; % Time span
M = eye(N); % Mass matrix
A_k = 5; % Amplitude of stiffness modulation
k_static = 500; % Baseline stiffness
k_wavenumber = 0.1*(2*pi); 
%k_wavenumber = 0;
k_angularfreq = 0.1*(2*pi);
%k_angularfreq = 0;
A_c = 0; % Amplitude of damping modulation
c_static = 0.1; % Baseline damping
c_wavenumber = 0*(2*pi);
c_angularfreq = 0*(2*pi);

B = zeros(N,1); % Forcing 
w_driving = 0;

%% Displacement and Filter
tic

a = GeneralizedForcedMSD(N, d, y0, ts, B, w_driving, M, A_k, k_static, k_wavenumber, k_angularfreq, A_c, c_static, c_wavenumber, c_angularfreq);
tiledlayout(3,2);


[t,y] = a.getStateVar();
U = y(:,1:N);   %Position data
V = y(:,N+1:end);   %Velocity data

% Single site displacements
nexttile([1 2])
plot(t,y(:,1),t,y(:,250));
axis on
title('Single Site Displacements')
xlabel('Time (s)') 
ylabel('Displacement (m)') 
legend({'Site 1','Site 250'},'Location','northeast')

% Single time displacements
nexttile([1 2])
plot(1:N, V(1000, :))
xlabel('Site No.')
ylabel(['Velocity at t = ' num2str(t(1000)) ])

% FFT
nexttile([1 2])
[f_1,Sig_1] = FFT_perso(t,y(:,1));
[f_250,Sig_250] = FFT_perso(t,y(:,250));
loglog(f_1,abs(Sig_1));
hold on
loglog(f_250,abs(Sig_250));

legend({'1st Mass FFT', '250th Mass FFT'},'Location','northeast');
title('Displacement FFT')
xlabel('frequency (Hz)') 
ylabel('amplitude') 



%% Eigenvalue change w/time
%nexttile([1 2])
figure
a_eigenvals = [];
for t_iter = 1:1:20
    [~,K] = a.getStiffness(t_iter);
    M_h = sqrt(M);
    K_tilde = M_h*K*M_h; % Mass normalized stiffness matrix
    [V,D] = eig(K_tilde);
    a_eigenvals = [a_eigenvals; diag(D)']
end

t_eigenvalues = 1:1:20;
for eigen_iter = 1:1:400
    plot(t_eigenvalues,a_eigenvals(:,eigen_iter));
    hold on
end
title('Eigenfrequency Evolution due to Stiffness Time Modulation')
xlabel('Time (s)') 
ylabel('Angular Frequency (rad/s)')      

%% Spatiotemporal Data
[SITES, TIMES] = meshgrid(1:N,t);    %Grid points for plots
U = y(:,1:N);   %Position data
V = y(:,N+1:end);   %Velocity data
figure
pcolor(SITES, TIMES, abs(V))
Ccolormap('Seahawks');
shading interp
caxis([0 0.05])
colorbar
xlabel('Site No.')
ylabel('Time')
title('|Velocity|')

%% Velocity vs Time base
% hold on
% 
% plot([0 N], t(2000)*[1 1], 'r--')
% plot(N/4*[1 1], [0 t(end)], 'm-.')
% 
% 
% figure
% subplot(2,1,1)
% plot(1:N, V(2000, :), 'ro-')
% xlabel('Site No.')
% ylabel(['Velocity at t = ' num2str(t(2000)) ])
% 
% subplot(2,1,2)
% plot(t, V(:,N/4), 'm-')
% xlabel('Time')
% ylabel(['Velocity of Site ', num2str(N/4)])

%% 2d FFT

f_cutoff = 2*sqrt(k_static/1); %characteristic frequency from dispersion relation 

dT = 2*pi/f_cutoff/100; %time step: 1/100*(period at cutoff frequency)

T_vec = dT*(0:4000);   %vector of output times

[FFT_V2, WaveNum2, Freq2] = FFT2_grid_v2(V, d, dT);

figure
pcolor(2*pi*WaveNum2, 2*pi*Freq2, fliplr(abs(FFT_V2)))

shading flat
Ccolormap('Seahawks')
% colorbar
% caxis([-60 0])
ylim([-8 8])
xlabel('\mu')
ylabel('\Omega')
%%
% %% Stiffness Angular Frequency Modulation Sweep
% k_mod_objects = [];
% for k_angularfreq = -0.2*(2*pi):0.1*(2*pi):0.2*(2*pi)
%         disp(k_angularfreq);
%         k_mod_objects = [k_mod_objects; GeneralizedForcedMSD(N, d, y0, ts, B, w_driving, M, A_k, k_static, k_wavenumber, k_angularfreq, A_c, c_static, c_wavenumber, c_angularfreq)];
% end
% %% 
% tiledlayout(2,1);
% nexttile
% for i = [1:5]
%     [t,y] = k_mod_objects(i).getStateVar();
%     plot(t,y(:,250));
%     hold on
% end
% legend({'-0.2 Hz','-0.1 Hz','0 Hz','0.1 Hz','0.2 Hz'},'Location','northeast');
% title('k angularfreq Sweep Displacements')
% xlabel('time (s)') 
% ylabel('displacement (m)') 
% 
% nexttile
% for i = [1:5]
%     [t,y] = k_mod_objects(i).getFilteredDisplacement(0,20,250);
%     [f,Sig] = FFT_perso(t,y); 
%     loglog(f,abs(Sig));
%     hold on
% end
% legend({'-0.2 Hz','-0.1 Hz','0 Hz','0.1 Hz','0.2 Hz'},'Location','northeast');
% title('k angularfreq Sweep FFT')
% xlabel('frequency (Hz)') 
% ylabel('amplitude') 
% %%
% close all
% clc
% 
% 
% N = 600; % Number of discrete masses/coordinates in the system
% d = .1; % Distance between masses
% 
% % Stiffness and damping modulation amplitude/rates
% A_k = 5;
% k_wavenumber = 0.1*(2*pi);
% k_angularfreq = 0.1*(2*pi);
% 
% A_c = 0;
% c_wavenumber = 0*(2*pi);
% c_angularfreq = 0*(2*pi);
% 
% generalizedForcedMSD(N, d, A_k, k_wavenumber, k_angularfreq, A_c, c_wavenumber, c_angularfreq)
% 
% %% State Variable Integration
% y0 = [0; zeros(N-1,1); 1; zeros(N-1,1)]; % Initial conditions, first half of vector is initial position and second half is initial velocity
% 
% ts = [0 20]; % Time domain (s)
% 
% %[t,y] = ode45('f',ts,y0);
%  options = odeset('AbsTol',1e-3,'RelTol',1e-3,'Stats','on');
% [t,y] = ode45('f',ts,y0,options);
% 
% 
% %% Plotting
% tiledlayout(6,2); % Requires R2019b or later
% 
% % Plot Initial k and c
% x = [1:1:N]*d;
% [c,C] = get_damping(0,k_wavenumber,k_angularfreq);
% [k,K] = get_stiffness(0,k_wavenumber,k_angularfreq);
% 
% nexttile
% scatter(x,k);
% title('Initial Stiffness Space-Dependence')
% xlabel('Distance (m)') 
% ylabel('Stiffness (N/m)') 
% 
% nexttile
% scatter(x,c);
% title('Initial Damping Space-Dependence')
% xlabel('Distance (m)') 
% ylabel('Damping (kg/s)') 
% 
% % Plot Displacement
% nexttile([1 2])
% plot(t,y(:,1),t,y(:,25));
% axis on
% title('Displacement vs. Time of Discrete Masses')
% xlabel('Time (s)') 
% ylabel('Displacement (m)') 
% legend({'Mass 1','Mass 25'},'Location','northeast')
% 
% 
% % Plot Velocity
% nexttile([1 2])
% plot(t,y(:,N+1),t,y(:,2*N));
% title('Velocity vs. Time of Discrete Masses')
% xlabel('Time (s)') 
% ylabel('Velocity (m/s)') 
% legend({'Mass 1','Mass 25'},'Location','northeast')
% axis on
% 
% %% FFT Analysis
% nexttile([1 2])
% 
% % Store displacement data
% m1_disp = y(:,1);
% m2_disp = y(:,25);
% 
% 
% % Create filtering window between x = filter_start and x = filter_end
% filter_start = 80;
% filter_end = 260;
% 
% w = [zeros(filter_start,1); window(@tukeywin,filter_end-filter_start); zeros(length(t)-filter_end,1)];
% m2_disp_filtered = y(:,25).*w;
% 
% plot(t,w,'r');
% hold on
% plot(t,m2_disp,'c');
% hold on
% plot(t,m2_disp_filtered,'b');
% 
% title('Windowing Function')
% xlabel('Time (s)') 
% ylabel('Displacement (m)') 
% 
% legend({'Windowing Function','Mass 25 Disp','Mass 25 Disp, Windowed'},'Location','northeast')
% 
% 
% nexttile([1 2])
% [f1,Sig1] = FFT_perso(t,m1_disp); % Perform FFT on 1st mass
% [f2,Sig2] = FFT_perso(t,m2_disp); % Perform FFT on 25th mass
% [f2_filt,Sig2_filt] = FFT_perso(t,m2_disp_filtered); % Perform FFT on 25th mass
% 
% 
% 
% plot(f1,abs(Sig1),'r')
% hold on
% plot(f2,abs(Sig2),'c')
% hold on
% plot(f2_filt,abs(Sig2_filt),'b')
% 
% legend({'Mass 1 (Raw)','Mass 25 (Raw)','Mass 25 (Windowed)'},'Location','northeast')
% xlim([0,4]);
% title('Displacement FFT')
% xlabel('Frequency (Hz)') 
% ylabel('Amplitude') 
% 
% 
% 
% %% Plot total system energy
%     % M = eye(N);
%     % E_displacement = 0.5*K.*(y(500,1:N).^2);
%     % E_velocity = 0.5*M.*(y(500,N+1:2*N));
% 
% [E_displacement, E_velocity, E_total] = getTotalEnergy(t,y,k);
% % 
% % nexttile
% % plot(t,E_displacement);
% % title('System Potential Energy')
% % xlabel('Time (s)') 
% % ylabel('U (J)') 
% % 
% % nexttile
% % plot(t,E_velocity)
% % title('System Kinetic Energy')
% % xlabel('Time (s)') 
% % ylabel('T (J)') 
% 
% nexttile([1 2])
% plot(t,E_total)
% title('System Total Energy')
% xlabel('Time (s)') 
% ylabel('Total Energy(J)') 
% 
% toc
% 
% %% Forward vs Reverse
% figure
% loglog(f_forward,sig_forward,'b')
% hold on
% loglog(f_reverse,sig_reverse,'r')
% title('Forward/Reverse FFT')
% xlabel('Frequency (Hz)') 
% ylabel('Amplitude') 
% 
% legend({'Mass 25 Forward','Mass 25 Reverse'},'Location','northeast')
% 
