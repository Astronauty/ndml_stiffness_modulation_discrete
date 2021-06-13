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
    end
    
    methods
        % Constructor
        function obj = GeneralizedForcedMSD(N, d, y0, ts, B, w_driving, M, A_k, k_static, k_wavenumber, k_angularfreq, A_c, c_static, c_wavenumber, c_angularfreq) 
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
            
            [obj.t,obj.y] = obj.integrateStateVar;

        end
        
        % Determines damping coefficient for each damper as a function of
        % space and time. c_output is a Nx1 vector of the dampers. C_output
        % is the corresponding NxN damping matrix.
        function [c_output,C_output] = getDamping(obj,t)
            C = eye(obj.N);
            c = zeros(obj.N,1);
            for i = 1:obj.N
                c(i) = obj.c_static+obj.A_c*sin(i*obj.d*obj.c_wavenumber-obj.c_angularfreq*t); % Damping coefficient of c = 0.1
            end

            % Create damping matrix C (tridiagonal) based on c
               C(obj.N,obj.N-1) = -c(obj.N-1); % Custom entries for last mass, since EQMS are different
               C(obj.N,obj.N) = c(obj.N);

             for row = 1:obj.N-1
                for col = 1:obj.N
                    if row == col % Diagonal elements
                        C(row,col) = c(row)+c(row+1);
                    elseif (row-col==-1) || (row-col==1)% Upper and lower diagonal elements (taking advantage of the fact the k index is tied to the column in both cases)
                        C(row,col) = -c(col);
                    end
                end
             end
     
             c_output = c;
             C_output = C;
        end
        
        function [k_output,K_output] = getStiffness(obj,t)
            K = eye(obj.N);
            
            % Stiffness vector k specifies initial spring rate at each discrete coordinate
            k = zeros(obj.N,1);
            for i = 1:obj.N
                % k(i) = mod(i-1, 5)+1;
                k(i) = obj.k_static+obj.A_k*sin(i*obj.d*obj.k_wavenumber-obj.k_angularfreq*t); % Baseline stiffness of 3, space modulation of 0.1 
            end

            % Create stiffness matrix K (tridiagonal) based on k
               K(obj.N,obj.N-1) = -k(obj.N-1); % Custom entries for last mass, since EQMS are different
               K(obj.N,obj.N) = k(obj.N);
                 for row = 1:obj.N-1
                    for col = 1:obj.N
                        if row == col % Diagonal elements
                            K(row,col) = k(row)+k(row+1);
                        elseif (row-col==-1) || (row-col==1)% Upper and lower diagonal elements (taking advantage of the fact the k index is tied to the column in both cases)
                            K(row,col) = -k(col);
                        end
                    end
                 end

               K_output = K;
               k_output = k;
        end
        
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
        
        function v=f(obj,t,x)
            % Damping Matrix C 
            [~,C] = obj.getDamping(t);

            % Stiffness Matrix K
            [~,K] = obj.getStiffness(t);

            % Return new state
            A1 = [zeros(obj.N) eye(obj.N); -obj.M\K -obj.M\C];
            %v = A1*x+[zeros(obj.N,1);f]*sin(obj.w_driving*t);
            f = obj.getForcing(t); 
            %f = zeros(obj.N,1);
            v = A1*x+[zeros(obj.N,1);f];
        end
        
        function [E_kinetic_out, E_potential_out, E_total_out] = getTotalEnergy(t,y,k)
           E_displacement = zeros(length(t),1);
           E_velocity = zeros(length(t),1);

           for i = 1:length(t)
               x = y(i,1:obj.N)';
               x_diff = x - [0;x(1:obj.N-1,1)];

               %test_displacementE = 0.5*K*(y(i,1:N)'.^2);
               E_displacement(i) = 0.5*sum(diag(k)*x_diff.^2);

               %E_displacement(i) = sum(0.5*K*(y(i,1:N)'.^2),'all');  
               E_velocity(i) = sum(0.5*obj.M*(y(i,obj.N+1:2*obj.N)'.^2),'all');
           end

           E_kinetic_out = E_displacement;
           E_potential_out = E_velocity;
           E_total_out = E_displacement+E_velocity;
        end
        
        function [T,Y] = integrateStateVar(obj)
            options = odeset('AbsTol',1e-3,'RelTol',1e-3,'Stats','on');
            [T,Y] = ode45(@(t,y) f(obj,t,y),obj.ts,obj.y0,options);
        end
        
        function [T,Y] = getStateVar(obj)
            [T,Y] = deal(obj.t, obj.y);
        end
        
        
        function [t,y] = getDisplacement(obj,mass_num)
            [t,y] = deal(obj.t, obj.y(:,mass_num));
        end
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
    end 
end


