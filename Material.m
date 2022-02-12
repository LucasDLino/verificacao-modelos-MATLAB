classdef Material < matlab.mixin.Copyable
    %Class Elastoplastic material
    %   Holds the properties of an elastoplastic material with Von Mises
    %   yield surface
    
    properties
        E
        nu
        K
        G
        stress
        strain
        elastic_strain
        plastic_strain
        effective_stress
        effective_strain
        equivalent_plastic_strain
        n_hard % number of sampling points for the piecewise linear hardening curve
        sampling_pairs
        H
        
        %Drucker-Prager model
        eta
        etabar
        xi
        
        %Mohr-Coulomb model
        phi
        psi
    end
    
    methods
        function obj = Material(young_modulus, poisson)
            %MATERIAL Constructor
            obj.E = young_modulus;
            obj.nu = poisson;
            obj.K = obj.E / (3. * (1. - 2. * obj.nu));
            obj.G = obj.E / (2. * (1. + obj.nu));
            
            %Defining other variables
            obj.stress = zeros(3);
            obj.strain = zeros(3);
            obj.elastic_strain = zeros(3);
            obj.plastic_strain = zeros(3);
            obj.effective_stress = 0;
            obj.effective_strain = 0;
            obj.equivalent_plastic_strain = 0;
            obj.sampling_pairs = [];
            obj.H = [];
            %Drucker-Prager
            obj.eta = 0;
            obj.etabar = 0;
            obj.xi = 0;
            %Mohr-Coulomb
            obj.phi = 0;
            obj.psi = 0;
        end
        
        function add_sampling_pair(obj, eps, sigmay)
            %Adds a sampling pair to the matrix
            obj.sampling_pairs = [obj.sampling_pairs; [eps, sigmay]];
            
            %Defining the hardening slope
            if size(obj.sampling_pairs, 1) >= 2
                hardening_slope = (obj.sampling_pairs(end, 2) - obj.sampling_pairs(end - 1, 2)) / (obj.sampling_pairs(end, 1) - obj.sampling_pairs(end - 1, 1));
                obj.H(end + 1) = hardening_slope;
            end
        end
        
        function delete(obj)
            %disp("Material deleted")
        end
    end
end

