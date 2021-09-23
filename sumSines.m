% Function to sum sinusoidal plane waves and generate a space/time mesh
% x - discrete positions to evaluate at 
% t - discete times to evaluate at
% w_x - wavenumber of plane wave (currently single value)
% w_t - vector which contains the frequencies you want to sum

function [X,T,outputWaveMesh] = sumSines(x,t,w_x,w_t_range)

    [X,T] = meshgrid(x,t);

    wavemesh = zeros(size(X));

    for w_t = w_t_range
        wavemesh = wavemesh + cos(w_x*X - w_t*T);
    end
    figure
    h = surf(X,T,wavemesh);
    %%h = plot(X,T,wavemesh);
    colormap winter;
    title("Stiffness Mesh");
    xlabel("Position (m)");
    ylabel("Time (s)");
    zlabel("Stiffness (N/m)");
    set(h, 'edgecolor','none')
    
    outputWaveMesh = wavemesh;
end

% function [x_out,t_out,k_out] = sumSines(x,t,A_k,k_static,w_x,w_t_range)
%     
%     k = zeros(x.length);
%     
%     i = 1:obj.N;
%     for w_t = w_t_range
%         k = k_static + A_k*(sin(w_x*X - w_t*T);
%         
%         k = obj.k_static+obj.A_k*(sin((i*obj.d)*obj.k_wavenumber - obj.k_angularfreq*t));
%     end
% 
%     h = surf(X,T,wavemesh);
%     %%h = plot(X,T,wavemesh);
%     colormap winter;
%     xlabel("Position (m)");
%     ylabel("Time (s)");
%     zlabel("Stiffness (N/m)");
%     set(h, 'edgecolor','none')
%     
%     outputWaveMesh = wavemesh;
% end