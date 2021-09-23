global k R L beta C m c C_array k_actual_array
R = 100;
L = 1E-6;
beta = 2e-8;
C = beta;
k = 1;
m = 1e-3;
c = 0.5;
C_array = [];
k_actual_array = [];


y0 = [0; 0; 0; 0]
[t,y] = ode15s(@(t,y) stateDerivative(y),[0 10], y0);

plot(t,y(:,1));
title("Displacement");
xlabel("Time (t)")
ylabel("Displacement (x)");

figure
plot(linspace(0,max(t),length(C_array)),(C_array));
title("Capacitance");
xlabel("Time (t)");
ylabel("Capacitance (F)")

figure
plot(linspace(0,max(t),length(C_array)),k_actual_array);
title("Effective Stiffness");
xlabel("Time (t)")
ylabel("Effective Stiffness (N/m)")

function z_dot = stateDerivative(z)
    global k R L beta C m c C_array k_actual_array
    x_equilibrium = 10e-3; % Temp solution for the fact that the state variable measures displacement from equil,but capacitance depends on absolute
    C = beta*(z(1)+x_equilibrium);
    C_array = [C_array C];
    v_in = 500;
    b = [0; (beta/(2*m))*z(4)^2; v_in/L; 0];

    k_actual_array = [k_actual_array (k+beta*z(4)^2/(z(1)+x_equilibrium))];
    %z_dot = Ax + b, but can't use Ax since nonlinear (capacitance term)
    A = [0 1 0 0;
         -(k/m) -(c/m) 0 0;
         0 0 -(R/L) -(1/L);
         0 0 1/C 0];

    z_dot = A*z + b;
end


