% Function for time and space dependent K matrix. Currently a plane wave
% solution

 function [c_output,C_output]=get_damping(t,c_wavenumber,c_angularfreq)
    global N d A_c
    
    C = eye(N);
    
    c = zeros(N,1);
    for i = 1:N
        c(i) = 0.1;
        %c(i) = 0+A_c*sin(i*d*c_wavenumber-c_angularfreq*t); % Damping coefficient of c = 0.1
    end
    
    % Create damping matrix C (tridiagonal) based on c
       C(N,N-1) = -c(N-1); % Custom entries for last mass, since EQMS are different
       C(N,N) = c(N);

     for row = 1:N-1
        for col = 1:N
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