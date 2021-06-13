% Function for time and space dependent K matrix. Currently a plane wave
% solution

 function [E_displacement_out, E_velocity_out, E_total_out] = getTotalEnergy(t,y,k)
    global N d A_k

   
   M = eye(N);
   E_total = zeros(length(t),1);
   E_displacement = zeros(length(t),1);
   E_velocity = zeros(length(t),1);
   
   for i = 1:length(t)
       x = y(i,1:N)';
       x_diff = x - [0;x(1:N-1,1)];

       %test_displacementE = 0.5*K*(y(i,1:N)'.^2);
       E_displacement(i) = 0.5*sum(diag(k)*x_diff.^2);
       
       %E_displacement(i) = sum(0.5*K*(y(i,1:N)'.^2),'all');  
       E_velocity(i) = sum(0.5*M*(y(i,N+1:2*N)'.^2),'all');
   end
  
   E_displacement_out = E_displacement;
   E_velocity_out = E_velocity;
   E_total_out = E_displacement+E_velocity;
 end