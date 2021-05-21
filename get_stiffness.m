% Function for time and space dependent K matrix. Currently a plane wave
% solution

 function [k_output,K_output]=get_stiffness(t,k_wavenumber,k_angularfreq)
    global N d A_k
    
    K = eye(N);
    
    % Stiffness vector k specifies initial spring rate at each discrete coordinate
    k = zeros(N,1);
    for i = 1:N
        % k(i) = mod(i-1, 5)+1;
        k(i) = 10+A_k*sin(i*d*k_wavenumber-k_angularfreq*t); % Baseline stiffness of 3, space modulation of 0.1 
    end
    
    
    
    % Create stiffness matrix K (tridiagonal) based on k
       K(N,N-1) = -k(N-1); % Custom entries for last mass, since EQMS are different
       K(N,N) = k(N);
         for row = 1:N-1
            for col = 1:N
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