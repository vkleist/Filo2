clear all
close all

param_set = {'Control';'Atg6_mutant'; 'Atg7_mutant';'Atg6_rescue'};

for i=1:length(param_set)
    mutant = char(param_set(i));
    switch mutant
        
        case 'Control'
            retract = 0;
            % Filopodia dynamics
            Parameters.(mutant).c1_sF = 1.82;%birth           
            Parameters.(mutant).c1_LF = 0.28;%birth
            Parameters.(mutant).c2_sF = 0.43;%death
            Parameters.(mutant).c2_LF = 0.07;%death
            %bulbous dynamics
            Parameters.(mutant).c3 = 0.0237;%
            Parameters.(mutant).B50 = 0.0776;%
            Parameters.(mutant).c4 = 1/120;%
            Parameters.(mutant).c5 = 0.0632;%
            % time dependent functions
            % a) dampening of filo dynamics
            Parameters.(mutant).p = [-2.9740e-17   3.3115e-13  -1.2896e-09   2.0637e-06  -1.4483e-03   1.0021e+00];
            % b) increase in propensity to form bulbous
            Parameters.(mutant).thalf = 1000;
            Parameters.(mutant).h = 1;
            
        case 'Atg6_mutant'
             retract = 0;%0.85;
            % Filopodia dynamics
            Parameters.(mutant).c1_sF = 1.19;%birth           
            Parameters.(mutant).c1_LF = 0.37;%birth
            Parameters.(mutant).c2_sF = 0.37;%death
            Parameters.(mutant).c2_LF = 0.05;%death
            %bulbous dynamics
            Parameters.(mutant).c3 = 0.0178;%0.2956;%
            Parameters.(mutant).B50 = 0.7157;%0.0584;%
            Parameters.(mutant).c4 = 1/120;%
            Parameters.(mutant).c5 = 0.0682;%0.2343;%
            % time dependent functions
            % a) dampening of filo dynamics
            Parameters.(mutant).p = Parameters.Control.p;
            % b) increase in propensity to form bulbous
            Parameters.(mutant).thalf = Parameters.Control.thalf;
            Parameters.(mutant).h = Parameters.Control.h;

        case 'Atg7_mutant'   
             retract = 0;%0.6;
            % Filopodia dynamics
            Parameters.(mutant).c1_sF = 1.48;%birth           
            Parameters.(mutant).c1_LF = 0.45;%birth
            Parameters.(mutant).c2_sF = 0.42;%death
            Parameters.(mutant).c2_LF = 0.05;%death
            %bulbous dynamics
            Parameters.(mutant).c3 = 0.0331;%0.1466;%0.0274;%
            Parameters.(mutant).B50 = 0.1618;%0.0489;%rB2S
            Parameters.(mutant).c4 = 1/120;%0.185;% %rB20
            Parameters.(mutant).c5 = 0.0751;%0.3549;%rB2S
            % time dependent functions
            % a) dampening of filo dynamics
            Parameters.(mutant).p = Parameters.Control.p;
            % b) increase in propensity to form bulbous
            Parameters.(mutant).thalf = Parameters.Control.thalf;
            Parameters.(mutant).h = Parameters.Control.h;       
            
            case 'Atg6_rescue'   
            retract = 0;%0.6;
            % Filopodia dynamics
            Parameters.(mutant).c1_sF = 2.1255;%birth           
            Parameters.(mutant).c1_LF = 0.1736;%birth
            Parameters.(mutant).c2_sF = 0.4545;%death
            Parameters.(mutant).c2_LF = 0.0769;%death
            %bulbous dynamics
            Parameters.(mutant).c3 = 0.0328;%0.1466;%0.0274;%
            Parameters.(mutant).B50 = 3.7331;%0.0489;%rB2S
            Parameters.(mutant).c4 = 0.0712;%0.185;% %rB20
            Parameters.(mutant).c5 = 0.0013;%0.3549;%rB2S
            % time dependent functions
            % a) dampening of filo dynamics
            Parameters.(mutant).p = Parameters.Control.p;
            % b) increase in propensity to form bulbous
            Parameters.(mutant).thalf = Parameters.Control.thalf;
            Parameters.(mutant).h = Parameters.Control.h;       
            
    end
end

try
    delete 'AllParameters.mat'
catch
end

save 'AllParameters.mat' Parameters 


