%% Model of Generalized 1D Discrete Nonreciprocal Vibration Transmission
% Daniel Nguyen
% 4/15/20
%% Unmodulated
clear all
tic

% Initialize shared variables
    n = 4; % Factorial samples

    N = 10;
    d = 1000E-6; % um
    y0 = [zeros(N,1); 0; zeros(N-1,1)]; % Initial conditions (currently 0 position, 0 velocity)
    ts = [0:0.00001:0.1]; % Time span
    M = eye(N); % Mass matrix
    
    k_static = linspace(1E6, 20E6, n); % Baseline stiffness
    k_base = linspace(1E6, 20E6, n);
    
    A_k = linspace(1E6,10E6, n); % Amplitude of stiffness modulation
    %k_wavenumber = 0.1*(2*pi); 
    k_wavenumber = 2*pi/(1E-3);
    k_angularfreq = linspace(1E2, 1E4, 20);
    
    A_c = 0; % Amplitude of damping modulation
    c_static = 0; % Baseline damping
    c_wavenumber = 0*(2*pi);
    c_angularfreq = 0*(2*pi);
    
    B = [10; zeros(N-1,1)]; % Forcing 
    w_driving = 100;
   
  factorial_matrix = fullfact([n n n 10]);
    a = GeneralizedForcedMSD(N, d, y0, ts, B, w_driving, M, 0 ...
        , k_static(factorial_matrix(1,2)), k_base((factorial_matrix(1,3))), k_wavenumber, k_angularfreq(factorial_matrix(1,4)), A_c, c_static, c_wavenumber, c_angularfreq);
    [t,y] = a.getStateVar();
    U = y(:,1:N);   %Position data
    V = y(:,N+1:end);   %Velocity data
    plot(t,y(:,5))
 
    %% Factorial Sample
    
    [row,col] = size(factorial_matrix);
    row*col
    pause
    candidates = [];
    candidates_indices = [];

    sampling_matrix = [A_k(factorial_matrix(:,1))' k_static(factorial_matrix(:,2))' k_base(factorial_matrix(:,3))' k_angularfreq(factorial_matrix(:,4))'];
    reciprocity_ratios = zeros(row,1);
    for i = 1:row
        forward = GeneralizedForcedMSD(N, d, y0, ts, B, w_driving, M, A_k(factorial_matrix(i,1))...
            , k_static(factorial_matrix(i,2)), k_base((factorial_matrix(i,3))), k_wavenumber, k_angularfreq(factorial_matrix(i,4)), A_c, c_static, c_wavenumber, c_angularfreq);
        
        [t_forward,y_forward] = forward.getStateVar();
        max_displacements_forward = max(y_forward(:,1:N));

        backward = GeneralizedForcedMSD(N, d, y0, ts, B, w_driving, M, A_k(factorial_matrix(i,1))...
            , k_static(factorial_matrix(i,2)), k_base((factorial_matrix(i,3))), k_wavenumber, -k_angularfreq(factorial_matrix(i,4)), A_c, c_static, c_wavenumber, c_angularfreq);
        [t_forward,y_backward] = backward.getStateVar();
        max_displacements_backward= max(y_backward(:,1:N));

        reciprocity_ratio = max_displacements_forward./max_displacements_backward;
        reciprocity_ratios(i) = max(reciprocity_ratio);
        
        if(abs(log(reciprocity_ratio)) > 2)
            candidates = [candidates; [A_k(factorial_matrix(i,1)) k_static(factorial_matrix(i,2)) k_base(factorial_matrix(i,3)) k_angularfreq(factorial_matrix(i,4))]]
            candidates_indices = [candidates_indices; i]
        end
        i
    end

    %%
    figure
    weights = reciprocity_ratios'*factorial_matrix
    
    scatter(sampling_matrix(:,1),reciprocity_ratios)

    set(gca,'yscale','log')
    xlabel("Stiffness Mod Amplitude (N/m*kg)")
    ylabel("Reciprocity Ratio")
   
    %%
    [coeff,score] = pca(candidates)
    figure
    vbls = {'Stiffness Mod Amplitude','Actuator-Actuator Stiffness','Onsite Stiffness',''};
    biplot(coeff(:,[1 2 3]),'Scores',score(:,[1 2 3]),'VarLabels',vbls)
