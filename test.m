%% Initialize default vars
tic
N = 20; % Number of masses
d = .005; % Distance between masses
y0 = [0; zeros(N-1,1); 1; zeros(N-1,1)]; % Initial conditions (currently 0 position, 0 velocity)
ts = [0 0.1]; % Time span
%M = 9E-6*eye(N); % Mass matrix
M = eye(N)
A_k = 420; % Amplitude of stiffness modulation
k_static = 21000; % Baseline stiffness
%k_wavenumber = 0.1*(2*pi); 
k_wavenumber = 0; 
%k_angularfreq = 0.1*(2*pi);
k_angularfreq = 0;
A_c = 0; % Amplitude of damping modulation
c_static = 0; % Baseline damping
c_wavenumber = 0*(2*pi);
c_angularfreq = 0*(2*pi);

B = zeros(N,1); % Forcing 
w_driving = 0;

%k_angularfreq_range = 5*pi:0.1:8*pi;
k_angularfreq_range = [5 25];
%%
a = GeneralizedForcedMSD(N, d, y0, ts, B, w_driving, M, A_k, k_static, k_wavenumber, k_angularfreq, A_c, c_static, c_wavenumber, c_angularfreq);

[t,y] = a.getStateVar();
U = y(:,1:N);   %Position data
V = y(:,N+1:end);   %Velocity data

[k,~] = a.getStiffness(10.1);
% figure
% plot(1:N,k);

%%
% % Single site displacements
% figure
% tiledlayout(1,2);
% nexttile([1 2])
% hold on
% plot(t,y(:,1));
% axis on
% title('Single Site Displacements');
% xlabel('Time (s)') 
% ylabel('Displacement (m)') 
% legend({'Site 1 Unmodulated'},'Location','northeast')

%%
figure
tiledlayout(1,2);
nexttile([1,2])

hold on
f_cutoff = 2*sqrt(k_static/M(1,1)); %cutoff frequency from dispersion relation 
dT = 2*pi/f_cutoff/100; %time step: 1/100*(period at cutoff frequency)
T_vec = dT*(0:4000);   %vector of output times
[FFT_V2, WaveNum2, Freq2] = FFT2_grid_v2(V, d, dT);

for n = 0
    pcolor((2*pi)*(WaveNum2)+n*k_wavenumber,(Freq2)- n*k_angularfreq, fliplr(abs(FFT_V2)));    
end




shading flat
Ccolormap('Seahawks')
% colorbar

xlim([-40 40])
ylim([-100 100])
xlabel('\mu')
ylabel('\Omega')
title('Unmodulated')

%%
% Variable definitions
syms wavenumber angularfreq f(w)

k = 500;
m = 10;
d = .1;

%Solve for dispersion relation
f(w) = sqrt(4*(k/m))*abs(sin(w*d/2));
angularfreq = f(wavenumber);

% Plot
fplot(f(wavenumber), [-pi/d pi/d],'r');
legend('Discrete Dispersion Relation','Analytical Dispersion Relation');
fplot(f(wavenumber-10)+2, [-pi/d pi/d]);
fplot(f(wavenumber+10)-2, [-pi/d pi/d]);
