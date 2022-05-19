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
        k_base
        k_wavenumber % Wavenumber stiffness modulation
        k_angularfreq % Angular freq stiffness modulation
        A_c % Amplitude of space damping modulation
        c_static 
        c_wavenumber
        c_angularfreq
        
        t
        y
        

    end
    
    methods
        % Constructor
        function obj = GeneralizedForcedMSD(N, d, y0, ts, B, w_driving, M, A_k, k_static, k_base, k_wavenumber, k_angularfreq, A_c, c_static, c_wavenumber, c_angularfreq) 
            obj.N = N;
            obj.d = d;
            obj.y0 = y0;
            obj.ts = ts;
            obj.B = B;
            obj.w_driving = w_driving;
            obj.M = M;
            obj.A_k = A_k;
            obj.k_static = k_static;
            obj.k_base = k_base
            obj.k_wavenumber = k_wavenumber;
            obj.k_angularfreq = k_angularfreq;
            obj.A_c = A_c;
            obj.c_static = c_static;
            obj.c_wavenumber = c_wavenumber;
            obj.c_angularfreq = c_angularfreq;
           
            
            [obj.t,obj.y] = obj.integrateStateVar;

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
             
%         function [k_output, K_output] = getStiffnessMultiple(obj,t)
%             k = zeros(1,obj.N);
%             i = 1:obj.N;
%             
%             for w_t = obj.k_angularfreq_range
%                 k = k + obj.A_k*sin((i*obj.d)*obj.k_wavenumber - w_t*t);
%             end    
%             %k = k/(length(obj.k_angularfreq_range)); % Normalize the amplitude based on how many frequencies summed
%             
%             k = k + obj.k_static;
%             K = diag(k)+diag([k(2:obj.N),0])+diag(-k(2:obj.N),1)+diag(-k(1:obj.N-1),-1);
% 
%             k_output = k;
%             K_output = K;
%             
%         end
        
        % Get stiffness vector k and stiffness matrix K at time t
        function [k_output,K_output] = getStiffness(obj,t)
            %K = zeros(obj.N);
            
            % Stiffness vector k specifies initial spring rate at each discrete coordinate
            i = 1:obj.N;
            k = zeros(1,obj.N);
            %k = obj.k_static+obj.A_k*(sin((i*obj.d)*obj.k_wavenumber - obj.k_angularfreq*t));
            
            for w_t = obj.k_angularfreq
                k = k + obj.A_k*sin((i*obj.d)*obj.k_wavenumber - w_t*t);
            end   

            % Create stiffness matrix K (tridiagonal) based on k
               k = k + obj.k_static;
               K = diag(k)+diag([k(1:obj.N)])+diag(-k(2:obj.N),1)+diag(-k(1:obj.N-1),-1);
               
               
               K_output = K + obj.k_base*eye(obj.N);
               k_output = k;
        end
        
        % Custom forcing function
        function f_out = getForcing(obj,t)
            % B not currently used for pulsed function implementation
            disp(t);
            f = zeros(obj.N,1);
            t_span = linspace(obj.ts(1),obj.ts(2),1000);
            
            t_w = linspace(0,2,1000)';
            w = window(@parzenwin,1000);
            
            if t > 0 && t < 2
                y = 100*sin(2*pi*4*t);
                y_filtered = y*w(1+round(500*t));
                
            else
                y_filtered = 0;
            end
            
            f(1) = y_filtered;

            f_out = f;
        end
        
        % State variable equation
        function v=f(obj,t,x)
            % Damping Matrix C 
            [~,C] = obj.getDamping(t);

            % Stiffness Matrix K
            %[~,K] = obj.getStiffness(t);
            [~,K] = obj.getStiffness(t);

            % Return new state
            A1 = [zeros(obj.N) eye(obj.N); -obj.M\K -obj.M\C];

%             v = A1*x+[zeros(obj.N,1);f]*sin(obj.w_driving*t);
            %f = obj.getForcing(t); 
%             f = zeros(obj.N,1);
            v = A1*x+[zeros(obj.N,1);obj.B];
        end

        % Get energies of the system 
        function [E_kinetic_out, E_potential_out, E_total_out] = getTotalEnergy(obj,t,y)
           E_potential = zeros(length(t),1);
           E_kinetic = zeros(length(t),1);
            
           [k,~] = obj.getStiffness(t);
           
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
            options = odeset('RelTol',1.e-3,'Stats','off');
%             [T,Y] = ode45(@(t,y) f(obj,t,y),obj.ts,obj.y0,options);
            [T,Y] = ode15s(@(t,y) f(obj,t,y),obj.ts,obj.y0,options);
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
                [k,~] = obj.getStiffness(i);
                x = obj.y(i,1:obj.N)'; % Displacements of each mass at time index i
                x_diff = (x - [0;x(1:obj.N-1,1)])'; % Find deflection of each spring based on differences of mass coordinates
                
                F(i,:) = k(i,:).*x_diff; % Spring forces, set a single row (site forces at given time)
                
            end
            F_out = F;
       
        end
        
    end 
end


