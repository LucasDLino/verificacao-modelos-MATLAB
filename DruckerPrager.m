classdef DruckerPrager < handle
    %DRUCKER-PRAGER MODEL class
    %   Responsible for updating the variables

    properties
        data_obj
    end

    methods
        function obj = DruckerPrager()
            %DRUCKER-PRAGER MODEL empty Constructor
        end

        function set_data(obj, input_data)
            %Sets up the data
            assert(isa(input_data,'Data'),'Incorrect data error:  input data is of class %s, not a Data object.', class(input_data));
            obj.data_obj = input_data;

            %Ajusting to analytical computation
            obj.data_obj.material_obj.sampling_pairs = [obj.data_obj.material_obj.sampling_pairs; [1e100, obj.data_obj.material_obj.sampling_pairs(end, 2)]];
            obj.data_obj.material_obj.H(end + 1) = 0;
        end

        function c = plfun(obj, cohesion)
            %Piecewise linear function defined by a set of pairs {ep, sy=f(ep)}

            sampling_pairs = obj.data_obj.material_obj.sampling_pairs;

            if obj.data_obj.material_obj.n_hard == 1
                c = sampling_pairs(1, 2);
                return
            end

            if cohesion < sampling_pairs(1, 1)
                c = sampling_pairs(1, 2);
                return
            end
            for i = 2:(obj.data_obj.material_obj.n_hard)
                if (cohesion < sampling_pairs(i, 1))
                    c = sampling_pairs(i - 1, 2) + ((cohesion - sampling_pairs(i - 1, 1)) * (sampling_pairs(i, 2) - sampling_pairs(i - 1, 2))) / (sampling_pairs(i, 1) - sampling_pairs(i - 1, 1));
                    return
                end
            end
            c = sampling_pairs(end, 2);
            return
        end

        function [is_plastic, is_fail] = computation(obj, strain)
            %Integration algorithm for the elastoplastic material with Drucker-Prager yield surface

            %Initialization of some algorithmic and internal variables
            dgama = 0; % incremental plastic multiplier
            tol = 10^(-12);
            is_plastic = false; % plastic yielding flag
            is_fail = false; % state update failure flag

            %Elastic predictor: Compute elastic trial state
            %----------------------------------------------------

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
            cohesion_trial = obj.plfun(eps_trial); % cohesion trial
            phi_trial = sqrt(varj2t_trial) + p_trial * obj.data_obj.material_obj.eta - obj.data_obj.material_obj.xi * cohesion_trial; % yield trial function
            res = phi_trial;

            if cohesion_trial ~= 0
                res = res/abs(cohesion_trial);
            end

            if (res > tol)
                %Plastic step: Apply return mapping to smooth portion of cone - computing analytically
                %----------------------------------------------------
                is_plastic = true;
                is_fail = true;

                eps = eps_trial;
                p = p_trial;
                s = s_trial;

                %Piecewise linear hardening
                for i = 1:obj.data_obj.material_obj.n_hard
                    dgama = (sqrt(varj2t_trial) + obj.data_obj.material_obj.eta * p_trial - obj.data_obj.material_obj.xi * (obj.data_obj.material_obj.sampling_pairs(i, 2) + obj.data_obj.material_obj.H(i) * ...
                        (eps_trial - obj.data_obj.material_obj.sampling_pairs(i, 1)))) / ...
                        (obj.data_obj.material_obj.G + obj.data_obj.material_obj.xi^2 * obj.data_obj.material_obj.H(i) + obj.data_obj.material_obj.eta * obj.data_obj.material_obj.etabar * obj.data_obj.material_obj.K);

                    eps = eps_trial + obj.data_obj.material_obj.xi * dgama;

                    if eps > obj.data_obj.material_obj.sampling_pairs(i, 1) && eps <= obj.data_obj.material_obj.sampling_pairs(i + 1, 1)
                        is_fail = false;
                        break;
                    end
                end

                p = p_trial - obj.data_obj.material_obj.K * obj.data_obj.material_obj.etabar * dgama;
                if varj2t_trial ~= 0
                    s = (1 - ((dgama*obj.data_obj.material_obj.G) / (sqrt(varj2t_trial)))) * s_trial;
                else
                    s = 0 * s_trial;
                end


                %Check validity of return to smooth portion
                if sqrt(varj2t_trial) - obj.data_obj.material_obj.G * dgama < 0
                    %Apply return mapping to apex portion of cone - computing analytically

                    %Set some variables
                    is_fail = true;
                    alpha = obj.data_obj.material_obj.xi / obj.data_obj.material_obj.etabar;
                    beta = obj.data_obj.material_obj.xi / obj.data_obj.material_obj.eta;

                    %Piecewise linear hardening
                    for i = 1:obj.data_obj.material_obj.n_hard
                        depv = (- beta * (obj.data_obj.material_obj.sampling_pairs(i, 2) + obj.data_obj.material_obj.H(i) * (eps_trial - obj.data_obj.material_obj.sampling_pairs(i, 1))) + p_trial) / ...
                            (alpha * beta * obj.data_obj.material_obj.H(i) + obj.data_obj.material_obj.K);

                        eps = eps_trial + alpha * depv;

                        if eps > obj.data_obj.material_obj.sampling_pairs(i, 1) && eps <= obj.data_obj.material_obj.sampling_pairs(i + 1, 1)
                            is_fail = false;
                            break;
                        end
                    end

                    p = p_trial - obj.data_obj.material_obj.K * depv;
                    s = 0 * s;
                end

                %Update state
                obj.data_obj.material_obj.strain = strain;
                obj.data_obj.material_obj.equivalent_plastic_strain = eps;
                obj.data_obj.material_obj.stress = s + (p * eye(3, 3));
                obj.data_obj.material_obj.effective_stress = sqrt(3. * 0.5 * trace(s * s));

                obj.data_obj.material_obj.elastic_strain = (s / (2. * obj.data_obj.material_obj.G)) + ((p * eye(3, 3)) / (3. * obj.data_obj.material_obj.K));
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

        function delete(~)
            %disp("Model deleted")
        end
    end
end

