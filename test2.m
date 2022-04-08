tic
N = 4;
%d = .005; % Distance between masses
d = .1;
y0 = [0; zeros(N-1,1); 10; zeros(N-1,1)]; % Initial conditions (currently 0 position, 0 velocity)
ts = [0 4]; % Time span
M = eye(N); % Mass matrix
k_static = 21000; % Baseline stiffness
A_k = 0.01*k_static; % Amplitude of stiffness modulation
k_wavenumber = 0.5*(2*pi); 
%k_wavenumber = 0;
k_angularfreq = 2*(2*pi);
%k_angularfreq = 0;
A_c = 0; % Amplitude of damping modulation
c_static = 20000; % Baseline damping
c_wavenumber = 0*(2*pi);
c_angularfreq = 0*(2*pi);

B = zeros(0,1); % Forcing 
w_driving = 0;

a = GeneralizedForcedMSD(N, d, y0, ts, B, w_driving, M, A_k, k_static, k_wavenumber, k_angularfreq, A_c, c_static, c_wavenumber, c_angularfreq);
[t,y] = a.getStateVar();

k_desired = [];
for i = 1:length(t)
    k_desired = [k_desired; a.getDesiredStiffness(y(i,1:N),t(i))];
end

displacement = y(:,1:N);   %Position data
velocity = y(:,N+1:end);   %Velocity data

differentialdisplacements = a.getDifferentialDisplacements();

actualStiffnesses = a.getDisplacementDependentStiffness(differentialdisplacements);

hold on
tiledlayout(2,1)

% Desired stiffness vs time plot
nexttile
plot(t,k_desired);
xlabel("Time (t)");
ylabel("Desired Stiffness (N/m)");

% Actual stiffness vs time plot
nexttile
plot(t,actualStiffnesses);
xlabel("Time (t)");
ylabel("Actual Stiffness (N/m)");
% Error in stiffness vs time plot

% nexttile
% %plot(t,displacement(:,1:2), 'DisplayName', 'mass1');
% plot(t,differentialdisplacements(:,1:4));
% 
% %plot(t,differential_displacement2, 'DisplayName', 'mass2');
% %plot(t,y(:,2), 'DisplayName', 'mass2');
% xlabel("Time (t)");
% ylabel("Differential Displacement (m)");
