%% Model of Generalized 1D Discrete Nonreciprocal Vibration Transmission
% Daniel Nguyen
% 4/11/20
%% Initialize default vars
tic
N = 500;
%d = .005; % Distance between masses
d = .1;
y0 = [0; zeros(N-1,1); 1; zeros(N-1,1)]; % Initial conditions (currently 0 position, 0 velocity)
ts = [0 2]; % Time span
M = eye(N); % Mass matrix
k_static = 21000; % Baseline stiffness
A_k = 0.01*k_static; % Amplitude of stiffness modulation
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

%% Displacement and Filter
tic

% Unmodulated solution
a = GeneralizedForcedMSD(N, d, y0, ts, B, w_driving, M, A_k, k_static, k_wavenumber, k_angularfreq, A_c, c_static, c_wavenumber, c_angularfreq);

[t,y] = a.getStateVar();
U = y(:,1:N);   %Position data
V = y(:,N+1:end);   %Velocity data

% Forward propagating wave
k_wavenumber = 2*(2*pi); 
k_angularfreq = (7:0.25:9)*(2*pi);
%k_angularfreq_range = [10*2*pi:pi/4:11*2*pi];
b = GeneralizedForcedMSD(N, d, y0, ts, B, w_driving, M, A_k, k_static, k_wavenumber, k_angularfreq, A_c, c_static, c_wavenumber, c_angularfreq);

[tb,yb] = b.getStateVar();
Ub = yb(:,1:N);
Vb = yb(:,N+1:end);

% Backward propagating wave
%k_angularfreq_range = -[10*2*pi:pi/4:11*2*pi];
k_angularfreq = -(7:0.25:9)*(2*pi);
c = GeneralizedForcedMSD(N, d, y0, ts, B, w_driving, M, A_k, k_static, k_wavenumber, k_angularfreq, A_c, c_static, c_wavenumber, c_angularfreq);


[tc,yc] = c.getStateVar();
Uc = yc(:,1:N);
Vc = yc(:,N+1:end);


%% Plotting
% Single site displacements
figure
tiledlayout(4,2);
nexttile([1 2])
hold on
plot(t,y(:,1),t,y(:,250),tb,yb(:,250),tc,yc(:,250));
axis on
title('Single Site Displacements');
xlabel('Time (s)') 
ylabel('Displacement (m)') 
legend({'Site 1 Unmodulated','Site 250 Unmodulated','Site 250 Forward Modulated','Site 250 Backward Modulated'},'Location','northeast')

% Single time displacements
nexttile([1 2])
plot(1:N, V(1, :),1:N, Vb(1, :),1:N, Vc(1, :))
xlabel('Site No.')
ylabel(['Velocity at t = ' num2str(t(1)) ])
legend({'Unmodulated','Forward Modulated', 'Backward Modulated'},'Location','northeast')


% FFT
nexttile([1 2])
hold on
[f_1,Sig_1] = FFT_perso(t,y(:,1));
[f_250,Sig_250] = FFT_perso(t,y(:,250));

[f_1b,Sig_1b] = FFT_perso(tb,yb(:,1));
[f_250b,Sig_250b] = FFT_perso(tb,yb(:,250));

[f_1c,Sig_1c] = FFT_perso(tc,yc(:,1));
[f_250c,Sig_250c] = FFT_perso(tc,yc(:,250));

hold on
loglog(f_1,abs(Sig_1));
loglog(f_250,abs(Sig_250));
loglog(f_1b,abs(Sig_1b));
loglog(f_250b,abs(Sig_250b));
loglog(f_1c,abs(Sig_1c));
loglog(f_250c,abs(Sig_250c));
hold off

legend({'1st Site Unmodulated', '250th Site Unmodulated','1st Site Forward Modulated','250th Site Forward Modulated','1st Site Backward Modulated','250th Site Backward Modulated' },'Location','northeast');
title('Displacement FFT')
xlim([0 30])
xlabel('frequency (Hz)') 
ylabel('amplitude') 

% Total Energy
nexttile([1 2])
[E_kinetic_a, E_potential_a, E_total_a] = a.getTotalEnergy(t,y);
[E_kinetic_b, E_potential_b, E_total_b] = b.getTotalEnergy(tb,yb);

hold on
plot(tb,E_total_b);
plot(tb,E_kinetic_b);
plot(tb,E_potential_b);
legend({'Total','Kinetic','Potential'},'Location','northeast');
title('System Total Energy');
xlabel('Time (s)');
ylabel('Total Energy(J)');


%% Spatiotemporal Data
figure
tiledlayout(3,1);
nexttile
[SITES, TIMES] = meshgrid(1:N,t);    %Grid points for plots
U = y(:,1:N);   %Position data
V = y(:,N+1:end);   %Velocity data
pcolor(SITES, TIMES, abs(V));
Ccolormap('Seahawks');
shading interp
caxis([0 0.05]);
colorbar
xlabel('Site No.')
ylabel('Time')
title('|Velocity|')


nexttile
[SITES, TIMES] = meshgrid(1:N,tb);    %Grid points for plots
Ub = yb(:,1:N);   %Position data
Vb = yb(:,N+1:end);   %Velocity data
pcolor(SITES, TIMES, abs(Vb));
Ccolormap('Seahawks');
shading interp
caxis([0 0.05]);
colorbar
xlabel('Site No.')
ylabel('Time')
title('|Velocity|')

nexttile
[SITES, TIMES] = meshgrid(1:N,tc);    %Grid points for plots
Uc = yc(:,1:N);   %Position data
Vc = yc(:,N+1:end);   %Velocity data
pcolor(SITES, TIMES, abs(Vc));
Ccolormap('Seahawks');
shading interp
caxis([0 0.05]);
colorbar
xlabel('Site No.')
ylabel('Time')
title('|Velocity|')

%% 2d FFT
figure
tiledlayout(2,2);

nexttile([1,2])
f_cutoff = 2*sqrt(k_static/M(1,1)); %characteristic frequency from dispersion relation 
dT = 2*pi/f_cutoff/100; %time step: 1/100*(period at cutoff frequency)
T_vec = dT*(0:4000);   %vector of output times
[FFT_V2, WaveNum2, Freq2] = FFT2_grid_v2(V, d, dT);

%surf(2*pi*WaveNum2, 2*pi*Freq2, fliplr(abs(FFT_V2)));
hold on
pcolor(2*pi*(WaveNum2),(Freq2), fliplr(abs(FFT_V2)));


syms wavenumber angularfreq f(w)


%Solve for dispersion relation
f(w) = sqrt(4*(k_static/M(1,1)))*abs(sin(w*d/2));
angularfreq = f(wavenumber);

% Plot analytical dispersion relation
fplot(f(wavenumber), [-pi/d pi/d],'r');
% fplot(f(wavenumber-abs(k_wavenumber/(2*pi)))+abs(k_angularfreq/(2*pi)), [-pi/d pi/d]);
% fplot(f(wavenumber+abs(k_wavenumber/(2*pi)))-abs(k_angularfreq/(2*pi)), [-pi/d pi/d]);
fplot(f(wavenumber-abs(k_wavenumber))+max(abs(k_angularfreq/(2*pi))), [-pi/d pi/d]);
fplot(f(wavenumber-abs(k_wavenumber))+min(abs(k_angularfreq/(2*pi))), [-pi/d pi/d]);
%fplot(f(wavenumber+k_wavenumber)-abs(10), [-pi/d pi/d]);

legend('Discrete Dispersion Relation','Analytical Dispersion Relation','Location','southeast');

% for w_t = k_angularfreq_range
%     surf(2*pi*(WaveNum2)+n*k_wavenumber, 2*pi*(Freq2)- w_t, fliplr(abs(FFT_V2)));
% end
view(2)

shading flat
Ccolormap('Seahawks')
% colorbar

xlim([-1.1*pi/d 1.1*pi/d])
ylim([-1.5*f_cutoff 1.5*f_cutoff])
xlabel('\mu')
ylabel('\Omega')
title('Unmodulated')


nexttile([1,1])
f_cutoff = 2*sqrt(k_static/M(1,1)); %characteristic frequency from dispersion relation 
dT = 2*pi/f_cutoff/100; %time step: 1/100*(period at cutoff frequency)
T_vec = dT*(0:4000);   %vector of output times
[backward_FFT_V2, backward_WaveNum2, backward_Freq2] = FFT2_grid_v2(Vc, d, dT);


pcolor(2*pi*backward_WaveNum2,-backward_Freq2, fliplr(abs(backward_FFT_V2)))
shading flat
Ccolormap('Seahawks')
% colorbar
caxis([1e-6 max(abs(FFT_V2),[],'all')])
xlim([-30 0])
ylim([-1.1*f_cutoff 1.1*f_cutoff])
xlabel('\mu')
ylabel('\Omega')
%title("Modulated Backward: " + "w_t=" + c.k_angularfreq + ", w_k=" + c.k_wavenumber);


nexttile([1,1])
f_cutoff = 2*sqrt(k_static/M(1,1)); %characteristic frequency from dispersion relation 
dT = 2*pi/f_cutoff/100; %time step: 1/100*(period at cutoff frequency)
T_vec = dT*(0:4000);   %vector of output times
[forward_FFT_V2, forward_WaveNum2, forward_Freq2] = FFT2_grid_v2(Vb, d, dT);


pcolor(2*pi*forward_WaveNum2, forward_Freq2, fliplr(abs(forward_FFT_V2)))
shading flat
Ccolormap('Seahawks')
% colorbar
caxis([1e-6 max(abs(forward_FFT_V2),[],'all')])
xlim([0 30])
ylim([-1.1*f_cutoff 1.1*f_cutoff])
xlabel('\mu')
ylabel('\Omega')
%title('Modulated Forward: ' + "w_t=" + b.k_angularfreq + ", w_k=" + b.k_wavenumber);

