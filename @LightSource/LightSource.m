classdef LightSource < handle
    %LigthSource Class for illuminant sources
    %
    
    properties
        ls_path
        ls
    end
    
    methods
        function obj = LightSource(p)
            %LightSource(p)
            %   Set the path to the light source and load the data
            
            % Set path
            obj.set_path(p);
            
            % Load data
            obj.load_data;
        end
        
        function obj = set_path(obj, p)
            %set_path
            %   Set path to light source
            obj.ls_path = p;
        end
        
        function obj = load_data(obj)
            %load_data
            % Load the data defined by ls_path
            
            load (obj.ls_path, 'spec');
            obj.ls = spec(1:10:401,2);
        end
        
    end
end

