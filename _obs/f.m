 function v=f(t,x)
    global N d k_wavenumber k_angularfreq c_wavenumber c_angularfreq
        % M = [2 0; 0 1]
        % C = [3 -0.5; -0.5 0.5]
        % K = [3 -1; -1 1]
        % B = [1; 1]
        % w = 1 % Driving frequency
 
     
    M = eye(N); % Mass matrix (kg)
    C = zeros(N); % Damping matrix (kg/s)
    B = zeros(N,1); % Forcing matrix (N) - using matrix rather than vector for future implementation of forcing modes if desired
    k = zeros(N,1);
    

    w = 2*pi; % Driving angular frequency (rad/s) - currently 1 Hz

    %% Mass Matrix M
    M = eye(N); % 1kg mass for every point


    %% Damping Matrix C 
    [c,C] = get_damping(t,c_wavenumber,c_angularfreq);

    %% Stiffness Matrix K
    [k,K] = get_stiffness(t,k_wavenumber,k_angularfreq);

    %% Forcing Matrix B
    B(1,1) = 0; % Only drive the first mass
    

    %% Return new state
%         A1 = [zeros(N) eye(N); -inv(M)*K -inv(M)*C];
%         f = inv(M)*B;
%         
        A1 = [zeros(N) eye(N); -M\K -M\C];
        f = M\B;
        v = A1*x+[zeros(N,1);f]*sin(w*t);
end
