% Variable definitions
syms wavenumber angularfreq f(w)

k = 21,000;
%k = 85.2; % N/m
%m = 0.00312; % kg
m = 9E-6;
d = .005; % m

%Solve for dispersion relation
f(w) = sqrt(4*(k/m))*abs(sin(w*d/2));
angularfreq = f(wavenumber);

% Plot
figure
hold
fplot(f(wavenumber), [-2*pi/d 2*pi/d]);
xlabel('\mu');
ylabel('\Omega');
fplot(f(wavenumber-20)+10, [-2*pi/d 2*pi/d]);
fplot(f(wavenumber+20)+10, [-2*pi/d 2*pi/d]);
% f(x,y) = x^2*y