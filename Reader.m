classdef Reader < handle
    %Reader class
    %   Responsible for reading the file and trasnfering variables to Data class
    
    properties
        file_name
        data_obj
        txt_in
    end
    
    methods
        function obj = Reader()
        end
        
        function set_fileName(obj, input_name)
            %Sets the name of the file
            obj.file_name = input_name;
            disp("File name: " + input_name);
        end
        
        function set_data(obj, input_data)
            %Sets the data
            assert(isa(input_data,'Data'),'Incorrect data error:  input data is of class %s, not a Data object.', class(input_data));
            obj.data_obj = input_data;
        end
        
        function read_elastic_properties(obj)
            %Reads the elastic properties
            tline = fgetl(obj.txt_in);
            splitline = str2double(split(tline));
            
            young = splitline(1);
            poisson = splitline(2);
            
            von_mises = Material(young, poisson);
            
            obj.data_obj.set_material(von_mises);
        end
        
        function read_vonMises_properties(obj)
            %Reads the elastoplastic hardening regions             
            pairs_number = str2double(fgetl(obj.txt_in));
            
            obj.data_obj.material_obj.n_hard = pairs_number;
            
            for i = 1:pairs_number
                pair = fgetl(obj.txt_in);
                splitpair = str2double(split(pair));
                
                eps = splitpair(1);
                sigmay = splitpair(2);
                
                obj.data_obj.material_obj.add_sampling_pair(eps, sigmay);
            end
        end
        
        function read_DruckerPrager_properties(obj)
            %Reads the elastoplastic hardening regions
            tline = fgetl(obj.txt_in);
            splitline = str2double(split(tline));
            
            pairs_number = splitline(1);
            
            obj.data_obj.material_obj.n_hard = pairs_number;
            
            obj.data_obj.material_obj.eta = splitline(2);
            obj.data_obj.material_obj.xi = splitline(3);
            obj.data_obj.material_obj.etabar = splitline(4);
            
            for i = 1:pairs_number
                pair = fgetl(obj.txt_in);
                splitpair = str2double(split(pair));
                
                eps = splitpair(1);
                c = splitpair(2);
                
                obj.data_obj.material_obj.add_sampling_pair(eps, c);
            end
        end
        
        function read_MohrCoulomb_properties(obj)
            %Reads the elastoplastic hardening regions
            tline = fgetl(obj.txt_in);
            splitline = str2double(split(tline));
            
            pairs_number = splitline(1);
            
            obj.data_obj.material_obj.n_hard = pairs_number;
            
            obj.data_obj.material_obj.phi = splitline(2);
            obj.data_obj.material_obj.psi = splitline(3);
            
            for i = 1:pairs_number
                pair = fgetl(obj.txt_in);
                splitpair = str2double(split(pair));
                
                eps = splitpair(1);
                c = splitpair(2);
                
                obj.data_obj.material_obj.add_sampling_pair(eps, c);
            end
        end
        
        function read_strain_history(obj)
            %Reads the strain history of the material           
            strain_number = str2double(fgetl(obj.txt_in));
            
            for i = 1:strain_number
                incsplit = str2double(split(fgetl(obj.txt_in)));
                
                strain = [
                    [incsplit(2), incsplit(5), incsplit(6)];
                    [incsplit(5), incsplit(3), incsplit(7)];
                    [incsplit(6), incsplit(7), incsplit(4)]
                    ];
                
                obj.data_obj.strains{end + 1} = strain;
            end
        end
        
        function read_all(obj)
            obj.txt_in = fopen(obj.file_name, "r");
                while ~feof(obj.txt_in)
                    tline = fgetl(obj.txt_in);
                    
                    if tline == "#ELASTIC.PROPERTIES"
                        obj.read_elastic_properties
                    end
                    
                    if tline == "#PLASTIC.VONMISES.PROPERTIES"
                        obj.read_vonMises_properties
                    end
                    
                    if tline == "#PLASTIC.DRUCKERPRAGER.PROPERTIES"
                        obj.read_DruckerPrager_properties
                    end
                    
                    if tline == "#PLASTIC.MOHRCOULOMB.PROPERTIES"
                        obj.read_MohrCoulomb_properties
                    end
                    
                    if tline == "#STRAIN.HISTORY"
                        obj.read_strain_history
                    end
                end
            fclose(obj.txt_in);
        end
        
        function delete(obj)
            %disp("Reader deleted")
        end
    end
end

