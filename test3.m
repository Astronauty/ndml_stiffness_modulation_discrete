tic
N = 4;
%d = .005; % Distance between masses
d = .005;
y0 = [0; zeros(N-1,1); 0; zeros(N-1,1)]; % Initial conditions (currently 0 position, 0 velocity)
ts = [0 10]; % Time span
M = 1*eye(N); % Mass matrix
k_static = 21000; % Baseline stiffness
A_k = 0.01*k_static; % Amplitude of stiffness modulation
%k_wavenumber = 0.02*(2*pi); 
k_wavenumber = 0;
%k_angularfreq = 0.5*(2*pi);
k_angularfreq = 0;
A_c = 0; % Amplitude of damping modulation
c_static = 2; % Baseline damping
c_wavenumber = 0*(2*pi);
c_angularfreq = 0*(2*pi);

B = 0;
w_drivings = 0;
hold on
for w_driving = w_drivings
    a = GeneralizedForcedMSD(N, d, y0, ts, B, w_driving, M, A_k, k_static, k_wavenumber, k_angularfreq, A_c, c_static, c_wavenumber, c_angularfreq);
    [t,y] = a.getStateVar();
    figure
    plot(t,y(:,2), 'DisplayName', 'mass1');
end
toc


