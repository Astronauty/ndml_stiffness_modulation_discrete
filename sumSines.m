% Function to sum sinusoidal plane waves and generate a space/time mesh
% x - discrete positions to evaluate at 
% t - discete times to evaluate at
% w_x - wavenumber of plane wave (currently single value)
% w_t - vector which contains the frequencies you want to sum
function outputWaveMesh = sumSines(x,t,w_x,w_t_range)

    [X,T] = meshgrid(x,t);

    wavemesh = zeros(size(X));

    % outputWave = zeros[
    % for w = k_cosineRange
    % sum(cos(k_wavenumber*x-w


    for w_t = w_t_range
        wavemesh = wavemesh + cos(w_x*X - w_t*T);
    end

    h = surf(X,T,wavemesh);
    colormap winter;
    xlabel("Position (m)");
    ylabel("Time (s)");
    zlabel("Stiffness (N/m)");
    set(h, 'edgecolor','none')
    outputWaveMesh = wavemesh;
    

end