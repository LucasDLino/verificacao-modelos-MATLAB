classdef Engine < handle
    %ENGINE class
    
    properties
        model_obj
        reader_obj
        data_obj
        output_obj
    end
    
    methods
        function obj = Engine()
            %ENGINE Constructor
            obj.model_obj = VonMises(); %Analytically
            %obj.model_obj = DruckerPrager(); %Analytically
            %obj.model_obj = MohrCoulomb(); %Analytically
            obj.reader_obj = Reader();
            obj.data_obj = Data();
            obj.output_obj = Output();
        end
        
        function start(obj, file_name)
            %Starts the simulation
            obj.reader_obj.set_fileName(file_name);            
            obj.reader_obj.set_data(obj.data_obj);
            obj.reader_obj.read_all;
            
            obj.model_obj.set_data(obj.data_obj);
            
            obj.output_obj.set_fileName(file_name);
            obj.output_obj.set_data(obj.data_obj);
            
            %obj.output_obj.store_state(false, false);
        end
        
        function solver(obj)
            %Solver: computation of the state variables history of the
            %material
            for k = obj.data_obj.strains
                [is_plastic, is_fail] = obj.model_obj.computation(k{1});
                obj.output_obj.store_state(is_plastic, is_fail);
            end
        end
        
        function finish(obj)
            %Finishes the simulation
            obj.output_obj.print_history;
            obj.output_obj.plot_history;
            %obj.output_obj.plot_axial;
            disp("Simulation ended successfully!");
        end
        
        function delete(obj)
            %disp("Engine deleted")
        end
    end
end

