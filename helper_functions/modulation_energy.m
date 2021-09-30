%% Measure Energy Required for Pump Wave
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

%%
% Forward propagating wave
k_wavenumber = 2*(2*pi); 
k_angularfreq = (7:0.25:9)*(2*pi);
b = GeneralizedForcedMSD(N, d, y0, ts, B, w_driving, M, A_k, k_static, k_wavenumber, k_angularfreq, A_c, c_static, c_wavenumber, c_angularfreq);

[tb,yb] = b.getStateVar();
Ub = yb(:,1:N);
Vb = yb(:,N+1:end);
%%
b.getModulationForce;

