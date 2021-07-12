%% Initialize default vars
tic
N = 500; % Number of masses
d = .1; % Distance between masses
y0 = [0; zeros(N-1,1); 1; zeros(N-1,1)]; % Initial conditions (currently 0 position, 0 velocity)
ts = [0 16]; % Time span
M = eye(N); % Mass matrix
A_k = 5; % Amplitude of stiffness modulation
k_static = 500; % Baseline stiffness
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

%%
a = GeneralizedForcedMSD(N, d, y0, ts, B, w_driving, M, A_k, k_static, k_wavenumber, k_angularfreq, A_c, c_static, c_wavenumber, c_angularfreq);
