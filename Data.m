classdef Data < handle
    %Data class
    %   Holds all data used in the simulation
    
    properties
        material_obj
        
        %Variable to store the state
        stresses
        effective_stresses
        effective_strains
        strains
        elastic_strains
        plastic_strains
        equivalent_plastic_strains
        is_plastic_step
        is_fail_step
    end
    
    methods
        function obj = Data()
            %Data - Construct an instance of this class  
            
            %Defining all variables
            obj.stresses = {};
            obj.effective_stresses = [];
            obj.effective_strains = [];
            obj.strains = {[0,0,0;0,0,0;0,0,0]};
            obj.elastic_strains = {};
            obj.plastic_strains = {};
            obj.equivalent_plastic_strains = [];
            obj.is_plastic_step = [];
            obj.is_fail_step = [];
        end
        
        function set_material(obj, input_material)
            %Sets up the material
            %   Veryfies if it is a material and sets the material
            assert(isa(input_material,'Material'),'Incorrect material error:  material is of class %s, not a Material object.', class(input_material));
            obj.material_obj = input_material;
        end
        
        function delete(obj)
            %disp("Data deleted")
        end
    end
end

