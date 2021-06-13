function map = Ccolormap(MAP,varargin)

%%% Custom colormap function. Each case below is a custom colormap, 
%%% and the otherwise block is for MATLAB built-ins.  

if ~ischar(MAP)
    error('Input MAP must be a string!')
end

N = 1;

if nargin > 1
    N = varargin{1};
end

switch MAP
    
    case 'Seahawks'
        map0 = ones(4*N,3);
        map0(:,1) = 2/3*ones(4*N,1);
        map0(:,2) = linspace(1, 0.9, 4*N);
        map0(:,3) = linspace(0.05, 0.4, 4*N);
        map0 = hsv2rgb(map0);

        map1 = ones(16*N,3);
        map1(:,1) = 2/3*ones(16*N,1);
        map1(:,2) = linspace(0.9,1,16*N);
        map1(:,3) = linspace(0.41, 0.90, 16*N);
        map1 = hsv2rgb(map1);

        map2 = ones(40*N,3);
        map2(:,3) = linspace(0.91, 0.95, 40*N);
        map2(:,2) = 1*map2(:,2);
        map2(:,1) = linspace(2/3, 1/3, 40*N);
        map2 = hsv2rgb(map2);

        map3 = ones(20*N,3);
        map3(:,1) = 1/3*ones(20*N,1);
        map3(:,2) = linspace(0.9,0.1,20*N);
        map3(:,3) = linspace(0.96,1,20*N);
        map3 = hsv2rgb(map3);

        map = [map0; map1; map2; map3];
        
    otherwise
        map = MAP;
        
end
        
colormap(map)

