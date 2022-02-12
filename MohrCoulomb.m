classdef MohrCoulomb < handle
    %MOHR-COULOMB MODEL class
    %   Responsible for updating the variables
    
    properties
        data_obj
    end
    
    methods
        function obj = MohrCoulomb()
            %MOHR-COULOMB MODEL empty Constructor
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
            %Integration algorithm for the elastoplastic material with Mohr-Coulomb yield surface
            
            %Initialization of some algorithmic and internal variables
            is_plastic = false; % plastic yielding flag
            is_fail = false; % state update failure flag
            
            %Getting material variables
            G = obj.data_obj.material_obj.G;
            K = obj.data_obj.material_obj.K;
            sampling_pairs = obj.data_obj.material_obj.sampling_pairs;
            H = obj.data_obj.material_obj.H;
            phi = obj.data_obj.material_obj.phi;
            psi = obj.data_obj.material_obj.psi;
            
            %Elastic predictor: Compute elastic trial state
           %----------------------------------------------------
            
            strain_increment = strain - obj.data_obj.material_obj.strain;
            
            %Defining the elastic trial state
            ee_trial = obj.data_obj.material_obj.elastic_strain + strain_increment; % trial elastic strain
            eps_trial = obj.data_obj.material_obj.equivalent_plastic_strain; % equivalent plastic trial strain
            eps = eps_trial;
            
            %Defining hydrostatic & volumetric stresses/strains
            eev_trial = trace(ee_trial); % elastric trial volumetric strain
            eed_trial = ee_trial - eev_trial * eye(3,3) / 3.; % elastic trial deviatoric strain
            
            p_trial = K * eev_trial; % hydrostatic stress
            s_trial = 2. * G * eed_trial; % deviatoric stress
            
            %Trial stress tensor
            stress_trial = s_trial + p_trial * eye(3,3);
            
            %Spectral decomposition of stress tensor
            [eigenvectors,eigenvalues] = eig(stress_trial); %Compute eigenvalues and eigenvectores
            [pstrs_trial, ind] = sort(diag(eigenvalues), 'descend'); %Sorted eigenvalues vector
            pstrs = pstrs_trial;
            %D_trial = eigenvalues(ind,ind); %Sorted eigenvalues diagonal matrix
            pdirs = eigenvectors(:,ind); %Sorted eigenvectors matrix
            
            %Check for plastic admissibility
            %----------------------------------------------------
            cohesion_trial = obj.plfun(eps_trial); % cohesion trial
            
            Phi_trial = pstrs_trial(1) - pstrs_trial(3) + (pstrs_trial(1) + pstrs_trial(3)) * sin(phi) - 2 * cohesion_trial * cos(phi); % yield trial function
            
            if (Phi_trial > 0)
                %Plastic step: Apply one-vector return mapping first (return to MAIN PLANE)
                %----------------------------------------------------
                is_plastic = true;
                
                dgamma = 0; % incremental plastic multiplier
                a = 4. * G * (1. + 1. / 3. * sin(phi) * sin(psi)) + 4. * K * sin(phi) * sin(psi);
                
                %Piecewise linear hardening - Explicit scheme
                is_fail = true;
                for i = 1:obj.data_obj.material_obj.n_hard
                    dgamma = ((sin(phi) + 1) * pstrs_trial(1) + (sin(phi) - 1) * pstrs_trial(3) - 2 * cos(phi) * (sampling_pairs(i,2) + (eps_trial - sampling_pairs(i,1)) * H(i))) /...
                        (a + 4 * H(i) * cos(phi) * cos(phi));
                    
                    eps = eps_trial + 2 * cos(phi) * dgamma;
                    
                    if eps >= sampling_pairs(i, 1) && eps <= sampling_pairs(i + 1, 1)
                        is_fail = false;
                        break;
                    end
                end
                
                pstrs(1) = pstrs_trial(1) - (2. * G * (1. + 1. / 3. * sin(psi)) + 2. * K * sin(psi)) * dgamma; % SIGMA1
                pstrs(2) = pstrs_trial(2) + (4. / 3. * G - 2. * K) * sin(psi) * dgamma; % SIGMA2
                pstrs(3) = pstrs_trial(3) + (2. * G * (1. - 1. / 3. * sin(psi)) - 2. * K * sin(psi)) * dgamma; % SIGMA3
                
                %TOL = 0;
                TOL = max(abs(pstrs)) * 10^-6;
                % Check validity of 1-vector return (check sextant of converged stress)
                if ~((pstrs(1) + TOL) >= pstrs(2) && (pstrs(2) + TOL) >= pstrs(3))
                    % Apply two-vector return mapping to appropriate EDGE
                    %----------------------------------------------------
                    %b = 0;
                    dgammaA = 0; dgammaB = 0;
                    
                    % Identify possible edge return: either right or left of main plane
                    edgeSide = pstrs_trial(1) * (1 - sin(psi)) - 2 * pstrs_trial(2) + pstrs_trial(3) * (1 + sin(psi));
                    
                    % Define sigmaA and sigmaB to simplify computation
                    sigmaA = pstrs_trial(1) - pstrs_trial(3) + (pstrs_trial(1) + pstrs_trial(3)) * sin(phi);
                    %sigmaB = 0;
                    
                    if (edgeSide > 0) % RIGHT Side
                        b = 2. * G * (1. + sin(phi) + sin(psi) - 1. / 3. * sin(phi) * sin(psi)) + 4. * K * sin(phi) * sin(psi);
                        sigmaB = pstrs_trial(1) - pstrs_trial(2) + (pstrs_trial(1) + pstrs_trial(2)) * sin(phi);
                        
                    else % LEFT Side
                        b = 2. * G * (1. - sin(phi) - sin(psi) - 1. / 3. * sin(phi) * sin(psi)) + 4. * K * sin(phi) * sin(psi);
                        sigmaB = pstrs_trial(2) - pstrs_trial(3) + (pstrs_trial(2) + pstrs_trial(3)) * sin(phi);
                    end
                    
                    % Explicit scheme for piecewise linear hardening
                    is_fail = true;
                    for i = 1:obj.data_obj.material_obj.n_hard
                        
                        dgammaA = (4. * H(i) * (sigmaA - sigmaB) * cos(phi) * cos(phi) - 2. * (a - b) * cos(phi) * (sampling_pairs(i,2) + (eps_trial - sampling_pairs(i,1)) * H(i)) + a * sigmaA - b * sigmaB) /...
                            ((a - b) * (8. * H(i) * cos(phi) * cos(phi) + a + b));
                        
                        dgammaB = (-4. * H(i) * (sigmaA - sigmaB) * cos(phi) * cos(phi) - 2. * (a - b) * cos(phi) * (sampling_pairs(i,2) + (eps_trial - sampling_pairs(i,1)) * H(i)) + a * sigmaB - b * sigmaA) /...
                            ((a - b) * (8. * H(i) * cos(phi) * cos(phi) + a + b));
                        
                        eps = eps_trial + 2 * cos(phi) * (dgammaA + dgammaB);
                        if (eps >= sampling_pairs(i,1) && eps <= sampling_pairs(i + 1,1))
                            is_fail = false;
                            break;
                        end
                    end
                    
                    if (edgeSide > 0)  % RIGHT Side
                        pstrs(1) = pstrs_trial(1) - (2. * G * (1. + 1. / 3. * sin(psi)) + 2. * K * sin(psi)) * (dgammaA + dgammaB); % SIGMA1
                        pstrs(2) = pstrs_trial(2) + (4. / 3. * G - 2. * K) * sin(psi) * dgammaA + (2. * G * (1. - 1. / 3. * sin(psi)) - 2. * K * sin(psi)) * dgammaB; % SIGMA2
                        pstrs(3) = pstrs_trial(3) + (2. * G * (1. - 1. / 3. * sin(psi)) - 2. * K * sin(psi)) * dgammaA + (4. / 3. * G - 2. * K) * sin(psi) * dgammaB; % SIGMA3
                        
                    else % LEFT Side
                        
                        pstrs(1) = pstrs_trial(1) - (2. * G * (1. + 1. / 3. * sin(psi)) + 2. * K * sin(psi)) * dgammaA + (4. / 3. * G - 2. * K) * sin(psi) * dgammaB; % SIGMA1
                        pstrs(2) = pstrs_trial(2) + (4. / 3. * G - 2. * K) * sin(psi) * dgammaA - (2. * G * (1. + 1. / 3. * sin(psi)) + 2. * K * sin(psi)) * dgammaB; % SIGMA2
                        pstrs(3) = pstrs_trial(3) + (2. * G * (1. - 1. / 3. * sin(psi)) - 2. * K * sin(psi)) * (dgammaA + dgammaB); % SIGMA3
                        
                    end
                end
                
                %TOL = 0;
                TOL = max(abs(pstrs)) * 10^-6;
                % Check validity of 2-vector return to edge
                if ~((pstrs(1) + TOL) >= pstrs(2) && (pstrs(2) + TOL) >= pstrs(3))
                    
                    % Apply multi-vector return mapping to APEX
                    %----------------------------------------------------
                    dEpv = 0;
                    alpha = cos(phi) / sin(psi);
                    
                    % Explicit scheme for piecewise linear hardening
                    is_fail = true;
                    for i = 1:obj.data_obj.material_obj.n_hard
                        
                        dEpv = (p_trial - cot(phi) * (sampling_pairs(i,2) + (eps_trial - sampling_pairs(i,1)) * H(i))) / (cot(phi) * alpha * H(i) + K);
                        eps = eps_trial + alpha * dEpv;
                        
                        if (eps >= sampling_pairs(i,1) && eps <= sampling_pairs(i + 1,1))
                            is_fail = false;
                            break;
                        end
                    end
                    
                    pstrs(:) = p_trial - K * dEpv;
                end
                
                %Update stress
                obj.data_obj.material_obj.stress = pdirs * diag(pstrs) * pdirs^(-1);
                
                p = trace(obj.data_obj.material_obj.stress)/3.;
                s = obj.data_obj.material_obj.stress - p * eye(3,3);
                
                %Update state for plastic step
                obj.data_obj.material_obj.strain = strain;
                obj.data_obj.material_obj.equivalent_plastic_strain = eps;
                
                obj.data_obj.material_obj.effective_stress = sqrt(3. * 0.5 * trace(s * s));

                obj.data_obj.material_obj.elastic_strain = (9. * K * obj.data_obj.material_obj.stress + trace(obj.data_obj.material_obj.stress) * eye(3,3) * (2. * G - 3. * K)) / (18. * G * K);
                %obj.data_obj.material_obj.elastic_strain = (s / (2. * G)) + ((p * eye(3, 3)) / (3. * K));
                
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

