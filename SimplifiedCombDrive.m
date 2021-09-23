classdef SimplifiedCombDrive
    %SIMPLIFIEDCOMBDRIVE Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
%         n % Number of teeth
%         b % Comb thickness
%         d % Comb spacing

        alpha % Proportionality between force and voltage: F = alpha*V^2
        V_in % Input voltage
        
        
    end
    
    methods
        function obj = SimplifiedCombDrive(alpha,V_in)
            obj.alpha = alpha;
            obj.V_in = V_in;
        end
        
        function F = getForce(obj)
            F = obj.alpha*obj.V_in^2;
        end
        
        
        function outputArg = method1(obj,inputArg)
            %METHOD1 Summary of this method goes here
            %   Detailed explanation goes here
            outputArg = obj.Property1 + inputArg;
        end
    end
end

