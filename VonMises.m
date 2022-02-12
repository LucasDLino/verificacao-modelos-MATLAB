classdef VonMises < handle
    %VON MISES MODEL class
    %   Responsible for updating the variables

    properties
        data_obj
    end

    methods
        function obj = VonMises()
            %VON MISES MODEL empty Constructor
        end

        function set_data(obj, input_data)
            %Sets up the data
            assert(isa(input_data,'Data'),'Incorrect data error:  input data is of class %s, not a Data object.', class(input_data));
            obj.data_obj = input_data;

            %Ajusting to analytical computation
            obj.data_obj.material_obj.sampling_pairs = [obj.data_obj.material_obj.sampling_pairs; [1e100, obj.data_obj.material_obj.sampling_pairs(end, 2)]];
            obj.data_obj.material_obj.H(end + 1) = 0;
        end

        function sigmay = plfun(obj, eps)
            %Piecewise linear function defined by a set of pairs {ep, sy=f(ep)}

            sampling_pairs = obj.data_obj.material_obj.sampling_pairs;

            if obj.data_obj.material_obj.n_hard == 1
                sigmay = sampling_pairs(1, 2);
                return
            end

            if eps < sampling_pairs(1, 1)
                sigmay = sampling_pairs(1, 2);
                return
            end
            for i = 2:(obj.data_obj.material_obj.n_hard)
                if (eps < sampling_pairs(i, 1))
                    sigmay = sampling_pairs(i - 1, 2) + ((eps - sampling_pairs(i - 1, 1)) * (sampling_pairs(i, 2) - sampling_pairs(i - 1, 2))) / (sampling_pairs(i, 1) - sampling_pairs(i - 1, 1));
                    return
                end
            end
            sigmay = sampling_pairs(end, 2);
            return
        end

        function [is_plastic, is_fail] = computation(obj, strain)
            %Integration algorithm for the elastoplastic material with Von Mises yield surface

            %Initialization of some algorithmic and internal variables
            is_plastic = false; % plastic yielding flag
            is_fail = false; % state update failure flag

            %Elastic predictor: Compute elastic trial state
            % ----------------------------------------------

            strain_increment = strain - obj.data_obj.material_obj.strain;

            %Defining the elastic trial state
            ee_trial = obj.data_obj.material_obj.elastic_strain + strain_increment; % trial elastic strain
            eps_trial = obj.data_obj.material_obj.equivalent_plastic_strain; % equivalent plastic trial strain

            %Defining hydrostatic & volumetric stresses/strains
            eev_trial = trace(ee_trial); % elastric trial volumetric strain
            eed_trial = ee_trial - eev_trial * eye(3,3) / 3.; % elastic trial deviatoric strain

            p_trial = obj.data_obj.material_obj.K * eev_trial; % hydrostatic stress
            s_trial = 2. * obj.data_obj.material_obj.G * eed_trial; % deviatoric stress

            %Check for plastic admissibility
            %----------------------------------------------------
            varj2t_trial = 0.5 * trace(s_trial * s_trial); % J2 invariant of the deviatoric stress tensor
            q_trial = sqrt(3 * varj2t_trial); % effective stress
            sigmay_trial = obj.plfun(eps_trial); % sigma_y trial
            phi_trial = q_trial - sigmay_trial; % yield trial function

            if (phi_trial > 0)
                %Plastic step: Apply return mapping
                %----------------------------------------------------
                dgama = 0; % incremental plastic multiplier
                is_plastic = true;
                is_fail = true;

                %Piecewise linear hardening
                for i = 1:obj.data_obj.material_obj.n_hard
                    dgama = (q_trial - (obj.data_obj.material_obj.sampling_pairs(i, 2) + obj.data_obj.material_obj.H(i) * (eps_trial - obj.data_obj.material_obj.sampling_pairs(i, 1)))) / (3. * obj.data_obj.material_obj.G + obj.data_obj.material_obj.H(i));

                    eps = eps_trial + dgama;

                    if eps >= obj.data_obj.material_obj.sampling_pairs(i, 1) && eps <= obj.data_obj.material_obj.sampling_pairs(i + 1, 1)
                        is_fail = false;
                        break;
                    end
                end

                %Update state
                obj.data_obj.material_obj.strain = strain;
                obj.data_obj.material_obj.equivalent_plastic_strain = eps_trial + dgama;
                s = (1 - ((dgama*3.*obj.data_obj.material_obj.G) / (q_trial))) * s_trial;
                obj.data_obj.material_obj.stress = s + (p_trial * eye(3, 3));
                obj.data_obj.material_obj.effective_stress = sqrt(3. * 0.5 * trace(s * s));
                obj.data_obj.material_obj.elastic_strain = (s / (2. * obj.data_obj.material_obj.G)) + (eev_trial * eye(3, 3) / 3.);
                obj.data_obj.material_obj.plastic_strain = obj.data_obj.material_obj.strain - obj.data_obj.material_obj.elastic_strain;

                ed = obj.data_obj.material_obj.strain - (trace(obj.data_obj.material_obj.strain)/3.) * eye(3,3);
                obj.data_obj.material_obj.effective_strain = sqrt((2. * trace(ed * ed)) / 3.);
            else
                %Elastic step: Update stress using linear elastic law
                %----------------------------------------------------
                obj.data_obj.material_obj.strain = strain;
                s = s_trial;
                obj.data_obj.material_obj.stress = s + (p_trial * eye(3, 3));
                obj.data_obj.material_obj.effective_stress = sqrt(3. * 0.5 * trace(s * s));
                obj.data_obj.material_obj.elastic_strain = ee_trial;
                obj.data_obj.material_obj.plastic_strain = obj.data_obj.material_obj.strain - obj.data_obj.material_obj.elastic_strain;
                obj.data_obj.material_obj.equivalent_plastic_strain = eps_trial;

                ed = obj.data_obj.material_obj.strain - (trace(obj.data_obj.material_obj.strain)/3.) * eye(3,3);
                obj.data_obj.material_obj.effective_strain = sqrt((2. * trace(ed * ed)) / 3.);
            end
        end

        function delete(obj)
            %disp("Model deleted")
        end
    end
end

