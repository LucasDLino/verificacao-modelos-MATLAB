classdef Output < handle
    %OUTPUT class
    %   Responsible for storing and printing the data
    
    properties
        data_obj
        file_name
        txt_out
    end
    
    methods
        function obj = Output()
            %OUTPUT empty constructor
        end
        
        function set_data(obj, input_data)
            %Sets the data
            assert(isa(input_data,'Data'),'Incorrect data error:  input data is of class %s, not a Data object.', class(input_data));
            obj.data_obj = input_data;
        end
        
        function set_fileName(obj, input_name)
            %Sets the name of the file
            obj.file_name = erase(input_name, ".txt") + "_output.txt";
            
            %Clear contents if file exists
            obj.txt_out = fopen(obj.file_name, "w");
            fclose(obj.txt_out);
            
%             if exist(obj.file_name, 'file') == 2
%                 delete(obj.file_name);
%             end 
        end
        
        function store_state(obj, is_plastic, is_fail)
            %Stores current state in data store variables
            obj.data_obj.stresses{end + 1} = obj.data_obj.material_obj.stress;
            obj.data_obj.effective_stresses(end + 1) = obj.data_obj.material_obj.effective_stress;
            obj.data_obj.effective_strains(end + 1) = obj.data_obj.material_obj.effective_strain;
            obj.data_obj.elastic_strains{end + 1} = obj.data_obj.material_obj.elastic_strain;
            obj.data_obj.plastic_strains{end + 1} = obj.data_obj.material_obj.plastic_strain;
            obj.data_obj.equivalent_plastic_strains(end + 1) = obj.data_obj.material_obj.equivalent_plastic_strain;
            obj.data_obj.is_plastic_step(end + 1) = is_plastic;
            obj.data_obj.is_fail_step(end + 1) = is_fail;
        end
        
        function print_history(obj)
            %Prints the history of the state variables of the simulation
            obj.txt_out = fopen(obj.file_name, "a");
            
            fprintf(obj.txt_out, "#STRESS.HISTORY\n");
            for stress = obj.data_obj.stresses
                state = [stress{1}(1,1), stress{1}(2,2), stress{1}(3,3), stress{1}(1,2), stress{1}(1,3), stress{1}(2,3)];
                fprintf(obj.txt_out, "%8.7f %8.7f %8.7f %8.7f %8.7f %8.7f\n", state);
            end
            
            fprintf(obj.txt_out, "\n");
            
            fprintf(obj.txt_out, "#STRAIN.HISTORY\n");
            for strain = obj.data_obj.strains
                state = [strain{1}(1,1), strain{1}(2,2), strain{1}(3,3), strain{1}(1,2), strain{1}(1,3), strain{1}(2,3)];
                fprintf(obj.txt_out, "%8.7f %8.7f %8.7f %8.7f %8.7f %8.7f\n", state);
            end
            
            fprintf(obj.txt_out, "\n");
            
            fprintf(obj.txt_out, "#EFFECTIVE.STRESS.HISTORY\n");
            for mises = obj.data_obj.effective_stresses
                fprintf(obj.txt_out, "%8.7f\n", mises);
            end
            
            fprintf(obj.txt_out, "\n");
            
            fprintf(obj.txt_out, "#EFFECTIVE.STRAIN.HISTORY\n");
            for ef_strain = obj.data_obj.effective_strains
                fprintf(obj.txt_out, "%8.7f\n", ef_strain);
            end
            
            fprintf(obj.txt_out, "\n");
            
            fprintf(obj.txt_out, "#ELASTIC.STRAIN.HISTORY\n");
            for elastic_strain = obj.data_obj.elastic_strains
                state = [elastic_strain{1}(1,1), elastic_strain{1}(2,2), elastic_strain{1}(3,3), elastic_strain{1}(1,2), elastic_strain{1}(1,3), elastic_strain{1}(2,3)];
                fprintf(obj.txt_out, "%8.7f %8.7f %8.7f %8.7f %8.7f %8.7f\n", state);
            end
            
            fprintf(obj.txt_out, "\n");
            
            fprintf(obj.txt_out, "#PLASTIC.STRAIN.HISTORY\n");
            for plastic_strain = obj.data_obj.plastic_strains
                state = [plastic_strain{1}(1,1), plastic_strain{1}(2,2), plastic_strain{1}(3,3), plastic_strain{1}(1,2), plastic_strain{1}(1,3), plastic_strain{1}(2,3)];
                fprintf(obj.txt_out, "%8.7f %8.7f %8.7f %8.7f %8.7f %8.7f\n", state);
            end
            
            fprintf(obj.txt_out, "\n");
            
            fprintf(obj.txt_out, "#EQUIVALENT.PLASTIC.STRAIN.HISTORY\n");
            for eps = obj.data_obj.equivalent_plastic_strains
                fprintf(obj.txt_out, "%8.7f\n", eps);
            end
            
            fprintf(obj.txt_out, "\n");
            
            fprintf(obj.txt_out, "#ISPLASTIC.ISFAIL.STEPS\n");
            logicalstr = {"false", "true"};
            for i = 1:length(obj.data_obj.is_plastic_step)
                state = [obj.data_obj.is_plastic_step(i) + 1, obj.data_obj.is_fail_step(i) + 1];
                fprintf(obj.txt_out, "%s %s\n", logicalstr{state(1)}, logicalstr{state(2)});
            end
            
            fclose(obj.txt_out);
        end
        
        function plot_history(obj)
            %Plots the loading history
           
            fig = figure("WindowState", "maximized");
            
            plot(obj.data_obj.effective_strains, obj.data_obj.effective_stresses,'-s',...
                'color', [0.6350 0.0780 0.1840],...
                'MarkerSize',5,...
                'MarkerFace',[0, 0.4470, 0.7410],...
                'MarkerEdgeColor',[0, 0.4470, 0.7410],...
                'LineWidth', 2);
            
            %title("Loading history", 'FontSize', 17);
            xlabel('Deformação equivalente', 'FontSize', 25);
            ylabel('Tensão equivalente (Pa)', 'FontSize', 25);
            
            labels = {};
            for i = 1:length(obj.data_obj.effective_strains)
                labels{end + 1} = i;
            end
            
            text(obj.data_obj.effective_strains, obj.data_obj.effective_stresses, cellfun(@num2str, labels,'un', 0),'FontSize', 12, 'VerticalAlignment','bottom','HorizontalAlignment','right');
            
            saveas(fig, erase(obj.file_name, ".txt") + ".png");
            saveas(fig, erase(obj.file_name, ".txt") + "_EPS.eps");
            pause(5);
        end
        
        function plot_axial(obj)
            stress_xx = [];
            strain_xx = [];
            labels = {};
            
            for i =1:length(obj.data_obj.stresses)
                stress_xx(end + 1) = obj.data_obj.stresses{i}(1,1);
                strain_xx(end + 1) = obj.data_obj.strains{i}(1,1);
                labels{end + 1} = int2str(i);
            end
            
            fig_axial = figure("WindowState", "maximized");
            
            plot(strain_xx, stress_xx,'-s',...
                'color', [0.6350 0.0780 0.1840],...
                'MarkerSize',5,...
                'MarkerFace',[0, 0.4470, 0.7410],...
                'MarkerEdgeColor',[0, 0.4470, 0.7410],...
                'LineWidth', 2);
            
            title("Loading history", 'FontSize', 17);
            xlabel('Strain_{xx}', 'FontSize', 12);
            ylabel('Stress_{xx} (Pa)', 'FontSize', 12);
            
            %text(strain_xx, stress_xx, labels, 'FontSize', 12, 'VerticalAlignment','bottom','HorizontalAlignment','right');
            
            saveas(fig_axial, erase(obj.file_name, ".txt") + "_axial.png");
        end
        
        function delete(obj)
            %disp("Output deleted")
        end
    end
end

