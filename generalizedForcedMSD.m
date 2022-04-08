classdef GeneralizedForcedMSD
    properties
        N % Number of discrete masses
        d % Distance between masses
        y0 % Initial conditions, first half of vector is initial position and second half is initial velocity
        ts % Time span
        B % Forcing Matrix
        M % Mass matrix
        w_driving % Driving frequency for forcing matrix
        A_k % Amplitude of space stiffness modulation
        k_static
        k_wavenumber % Wavenumber stiffness modulation
        k_angularfreq % Angular freq stiffness modulation
        A_c % Amplitude of space damping modulation
        c_static 
        c_wavenumber
        c_angularfreq
        
        t
        y
        
        error
        error_derivative
        error_integral

    end
    
    methods
        % Constructor
        function obj = GeneralizedForcedMSD(N, d, y0, ts, B, w_driving, M, A_k, k_static, k_wavenumber, k_angularfreq, A_c, c_static, c_wavenumber, c_angularfreq); 
            obj.N = N;
            obj.d = d;
            obj.y0 = y0;
            obj.ts = ts;
            obj.B = B;
            obj.w_driving = w_driving;
            obj.M = M;
            obj.A_k = A_k;
            obj.k_static = k_static;
            obj.k_wavenumber = k_wavenumber;
            obj.k_angularfreq = k_angularfreq;
            obj.A_c = A_c;
            obj.c_static = c_static;
            obj.c_wavenumber = c_wavenumber;
            obj.c_angularfreq = c_angularfreq;
           
            % PID Terms
            obj.error = [];
            obj.error_integral = [];
            obj.error_derivative = [];
            
            [obj.t,obj.y] = obj.integrateStateVar;

        end
        
        function [k_output] = getDisplacementDependentStiffness(obj, y)
            %k_output = exp(500*y)' + obj.k_static;
            k_output = 2000*y' + obj.k_static;
        end
        
        function [displacement] = getStiffnessDependentDisplacement(obj, k)
            %k_output = exp(500*y)' + obj.k_static;
            displacement = [(k - obj.k_static)./2000.0]';
        end
        
        function differential_displacements = getDifferentialDisplacements(obj)
            displacements = obj.y(:,1:obj.N);
            [rows, cols] = size(displacements);
            differential_displacements = displacements - [zeros(rows, 1) displacements(:,1:(cols-1))];
        end
        
        % Determines damping coefficient for each damper as a function of
        % space and time. c_output is a Nx1 vector of the dampers. C_output
        % is the corresponding NxN damping matrix.
        function [c_output,C_output] = getDamping(obj,t)
            i = 1:obj.N;
            c = obj.c_static+obj.A_c*sin(i*obj.d*obj.c_wavenumber-obj.c_angularfreq*t);

            % Create damping matrix C (tridiagonal) based on c
             C = diag(c)+diag([c(2:obj.N),0])+diag(-c(2:obj.N),1)+diag(-c(1:obj.N-1),-1);

             c_output = c;
             C_output = C;
        end
             
        function [k_output,K_output] = getDesiredStiffness(obj,y,t)
            % Stiffness vector k specifies initial spring rate at each discrete coordinate
            i = 1:obj.N;
            k = zeros(1,obj.N);

            % Add up each of the angular frequency components of modulation
            % ASSUMES SMALL ENOUGH DISPLACEMENT THAT IT DOESN'T CONTRIBUTE
            % TO DISPLACEMENT DEPENDENT STIFFNESS CHANGE
            for w_t = obj.k_angularfreq
                k = k + obj.A_k*sin((i*obj.d)*obj.k_wavenumber - w_t*t);
            end   
            
            k = k+obj.k_static;
            % Create stiffness matrix K (tridiagonal) based on k
               K = diag(k)+diag([k(2:obj.N),0])+diag(-k(2:obj.N),1)+diag(-k(1:obj.N-1),-1);
               
               K_output = K;
               k_output = k;
        end

        % Get stiffness vector k and stiffness matrix K at time t
        function [k_output,K_output] = getActualStiffness(obj,y,t)
            % Stiffness vector k specifies initial spring rate at each discrete coordinate
            i = 1:obj.N;
            k = zeros(1,obj.N);

            % Add up each of the angular frequency components of modulation
            for w_t = obj.k_angularfreq
                k = k + obj.A_k*sin((i*obj.d)*obj.k_wavenumber - w_t*t);
            end   
            
            k = k + obj.getDisplacementDependentStiffness(y(1:obj.N)); % Add the baseline stiffness
            % k = k+obj.k_static;
            % Create stiffness matrix K (tridiagonal) based on k
               K = diag(k)+diag([k(2:obj.N),0])+diag(-k(2:obj.N),1)+diag(-k(1:obj.N-1),-1);
               
               K_output = K;
               k_output = k;
        end
        
        % Get energies of the system 
        function [E_kinetic_out, E_potential_out, E_total_out] = getTotalEnergy(obj,t,y)
           E_potential = zeros(length(t),1);
           E_kinetic = zeros(length(t),1);
            
           [k,~] = obj.getActualStiffness(t);
           
           for i = 1:length(t)
               x = y(i,1:obj.N)'; % Displacements of each mass at time index i
               x_diff = x - [0;x(1:obj.N-1,1)]; % Find deflection of each spring based on differences of mass coordinates

               %test_displacementE = 0.5*K*(y(i,1:N)'.^2);
               E_potential(i) = 0.5*sum(k(i,:)*x_diff.^2); % Spring potential energy

               %E_displacement(i) = sum(0.5*K*(y(i,1:N)'.^2),'all');  
               E_kinetic(i) = sum(0.5*obj.M*(y(i,obj.N+1:2*obj.N)'.^2),'all'); % Kinetic energy
           end

           E_kinetic_out = E_potential;
           E_potential_out = E_kinetic;
           E_total_out = E_potential+E_kinetic;
        end
        
        function [T,Y] = integrateStateVar(obj)
            %options = odeset('RelTol',1.e-4,'Stats','on');
            options = odeset('RelTol', 1e-3, 'MaxStep',0.1);
           
            [T,Y] = ode15s(@(t,y) f(obj,t,y),obj.ts,obj.y0,options);
        end
        
                % State variable equation
        function v=f(obj,t,y)
            % Damping Matrix C 
            [~,C] = obj.getDamping(t);

            % Stiffness Matrix K
            [k,K] = obj.getActualStiffness(y,t);
            
            % Solve for PID forcing
            [k_desired,~] = obj.getDesiredStiffness(y,t);
            desiredDifferentialDisplacement = obj.getStiffnessDependentDisplacement(k_desired);
            
            desiredDisplacement = zeros(obj.N,1);
            
            desiredDisplacement(1) = desiredDifferentialDisplacement(1);
            for i = 2:length(desiredDifferentialDisplacement)
                desiredDisplacement(i) = desiredDisplacement(i-1) + desiredDifferentialDisplacement(i);
            end
                
            
            PID_force = obj.getPIDForce(desiredDisplacement, y(1:obj.N));
            % Return new state
            A1 = [zeros(obj.N) eye(obj.N); -obj.M\K -obj.M\C];
            
            f = [1; 0; 0; 0];
%            v = A1*y+[zeros(obj.N,1);f*sin(obj.w_driving*t)];

            %f = zeros(obj.N,1);
            %f = [1000; 1000; -500; 1000]*sin(8*t);
            f = PID_force;
            v = A1*y+[zeros(obj.N,1);f];
        end

        function f = getPIDForce(obj, desiredDisplacement, actualDisplacement)
            % Rows are the displacement of each mass, columns are each
            % time
            error = desiredDisplacement - actualDisplacement;
            obj.error = [obj.error error];
            Kp = 500000;
            f = Kp * error;
        end
        

        

        % Get the time-displacement data 
        function [T,Y] = getStateVar(obj)
            [T,Y] = deal(obj.t, obj.y);
        end
        
        % Get the displacement of a particular mass/site
        function [t,y] = getDisplacement(obj,mass_num)
            [t,y] = deal(obj.t, obj.y(:,mass_num));
        end
        
        % Get displacement after a windowing function has been applied
        function [t,y_filtered] = getFilteredDisplacement(obj,t_start,t_end,mass_num)
            i = 1;
            while obj.t(i) < t_start
                i = i+1;
            end
            start_index = i;
            
            i = 1;
            while obj.t(i) < t_end
                i = i+1;
            end
            end_index = i;
            
            w = [zeros(start_index,1); window(@tukeywin,end_index-start_index); zeros(length(obj.t)-end_index,1)];
            
            t = obj.t;
            y_filtered = obj.y(:,mass_num).*w;
            
        end
        
        % Get the modulation forces across all of the sites at a given time
        function [F_out] = getModulationForce(obj)
            F = zeros(length(obj.t),obj.N);
            
            for i = 1:length(obj.t)
                [k,~] = obj.getActualStiffness(i);
                x = obj.y(i,1:obj.N)'; % Displacements of each mass at time index i
                x_diff = (x - [0;x(1:obj.N-1,1)])'; % Find deflection of each spring based on differences of mass coordinates
                
                F(i,:) = k(i,:).*x_diff; % Spring forces, set a single row (site forces at given time)
                
            end
            F_out = F;
       
        end
        
    end 
end


