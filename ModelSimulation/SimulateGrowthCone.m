%growthConeSim
close all
clear all

%Gillespie ensemble setup
N_gil = 50;%number of simulations
T = 3600;%final time (60 hours = 3600 minutes = 216000 seconds)
AutophagyType = 'Atg6_rescue';%'Control';'Atg6_mutant'; 'Atg7_mutant'; 'Atg6_rescue'

%output file
file = strcat('EnsembleData/DataSimple_',AutophagyType);

% Options for saving ensemble data
time_disc = 1;%take a snapshot every minute
time_points = ceil(T/time_disc);
ens_Data  = zeros(N_gil,time_points,5); %ensemble data filopodia

[c1_sF,c1_LF,c2_sF,c2_LF,c3, B50, c4, c5, p, thalf, h] = getParameters(AutophagyType);% c = c1, c2, growth_rate_F, d_F, k_sF, d_sF, k_B, k_B2, d_B, B_50


Initial_sFilo = c1_sF/c2_sF;
Initial_LFilo = c1_LF/c2_LF;

% Stoichiometric matrix     
%   r1a r1b r2a r2b r3 r4a r4b  rsBLB   r5    
S = [1  0   -1  0   0   0   0     0      0;%sF
     0  1   0   -1  0   0   0     0      0;%LF
     0  0   0   0   1   -1  0     -1      0;%sB
     0  0   0   0   0   0   -1    1      -1;%LB
     0  0   0   0   0   0   0     0      1];%S

%start of Gillespie run
for i=1:N_gil
tic
t = 0; %start time

%Initialize filopodia/stab filo/bilbous/synapse
sF = poissrndMat(Initial_sFilo); %short lived filopodia
LF = poissrndMat(Initial_LFilo); %long lived folopodia
%actual state
X = [sF;LF;0;0;0];
sBul_LifeTimes = nan(100,1);

%setup storage
counter2 = 0;
SizeStoreVector = 5000;
Time = nan(1,SizeStoreVector);
XStore = nan(length(X),SizeStoreVector);

% Update state
counter2= counter2+1;
Time(counter2) = 0;
XStore(:,counter2) = X;

counter = 0;
%start of growth cone simulation
while t < T
    
    %compute propensity and their sum
    Vec = sBul_LifeTimes(~isnan(sBul_LifeTimes));
    if isempty(Vec)
        Vec = -inf;
    end
    a = getPropensities(t,X,c1_sF,c1_LF,c2_sF,c2_LF,c3, B50, c4, c5, thalf, h,p);
    a0 = sum(a);
    % 1. time update
    r1 = rand;
    tau = (1/a0)*log(1/r1);%time step
    t = t+tau; %time update
    % 2. sample reaction
    r2 = rand;
    j = find(r2 <= cumsum(a)/a0,1);
    
    X_before = X;
    if X(3) < 0
        disp('before update')
        return
    end
  
    % 3. state update
    X = X + S(:,j);
    
    if X(3) < 0
        return
    end

    % reaction R4a, sB died
    if j == 6 
        idx = find(~isnan(sBul_LifeTimes),1,'last');
        sBul_LifeTimes(idx) = nan;%update lifetimes
    end
    
    %Update Bulbous Lifetimes
    sBul_LifeTimes = sBul_LifeTimes+tau;
    
    %emergence of a new bulbous; reaction R3, initialize its lifetime
    if j == 5
        idx = find(isnan(sBul_LifeTimes),1);
        sBul_LifeTimes(idx) = 0;
    end
    
    %check whether sB -> LB
    if j == 8
        Vec = sBul_LifeTimes(~isnan(sBul_LifeTimes));
        r3 = rand;
        idx2 = find(~isnan(sBul_LifeTimes));
        z = find(r3 <= cumsum(Vec)/sum(Vec),1);
        sBul_LifeTimes(idx2(z)) = nan;
    end

    
    %time data for ensemble matrices
    t_now = ceil(t);
    t_before = ceil(Time(counter2));
    t_passed = t_before:time_disc:t_now-1;
    
    if t_now > T
        break;
    end
    
    %Matrices updates
    counter2 = counter2+1;
    if counter2 > SizeStoreVector
        Time = [Time,nan(1,5000)];
        XStore= [ XStore,nan(5,5000)];
        SizeStoreVector = SizeStoreVector+5000;
    end
    Time(counter2) = t;
    XStore(:,counter2) = X;
    
    %update ensemble matrices
    ens_Data(i,t_now,1) = X(1);
    ens_Data(i,t_now,2) = X(2);
    ens_Data(i,t_now,3) = X(3);
    ens_Data(i,t_now,4) = X(4);
    ens_Data(i,t_now,5) = X(5);

    
    if t_before <= 1
        continue;
    end
    
    %if t_passed isn't 0, fill in data for the time between
    aux = ens_Data(i,t_before,1);
    ens_Data(i,t_passed,1) = aux;
    %
    aux = ens_Data(i,t_before,2);
    ens_Data(i,t_passed,2) = aux;
    %
    aux = ens_Data(i,t_before,3);
    ens_Data(i,t_passed,3) = aux;
    %
    aux = ens_Data(i,t_before,4);
    ens_Data(i,t_passed,4) = aux;
    %
    aux = ens_Data(i,t_before,5);
    ens_Data(i,t_passed,5) = aux;
    

    
end

toc

%fill in last bits of ensemble matrices
t_before = ceil(Time(counter2));
t_passed = t_before:time_disc:T;

aux = ens_Data(i,t_before,1);
ens_Data(i,t_passed,1) = aux;
%
aux = ens_Data(i,t_before,2);
ens_Data(i,t_passed,2) = aux;
%
aux = ens_Data(i,t_before,3);
ens_Data(i,t_passed,3) = aux;
%
aux = ens_Data(i,t_before,4);
ens_Data(i,t_passed,4) = aux;
%
aux = ens_Data(i,t_before,5);
ens_Data(i,t_passed,5) = aux;
end

save(file,'ens_Data');
%Trimming
Time = Time(1:counter2);
XStore = XStore(:,1:counter2);

PlotSimulation(Time,XStore,time_points,ens_Data,AutophagyType);


function a = getPropensities(t,X,c1_sF,c1_LF,c2_sF,c2_LF,c3, B50, c4, c5, thalf, h,p)

    a = zeros(1,9);
    a(1)  = c1_sF*time_dampening(t,p); %birth of short lived filo
    a(2)  = c1_LF*time_dampening(t,p); %birth of long lived filo
    a(3)  = c2_sF*X(1); %death of short lived filopodium
    a(4)  = c2_LF*X(2); %death of long lived filopodium
    a(6)  = c4*(X(3)); % short-bulbous death
    a(7)  = 0; % long-bulbous death
    a(9)  = 1/133*X(4);%; % bulbous to synapse
    
    a(5)  = c3*(X(1)+X(2))*time_enhancing(t,thalf,h); % Filopodium-to-bulbous transition
    a(8)  = X(3)*c5*feedback(B50,X(4));%short-to-long bulb

end

%Function that influences the enhancement of bulbous birth
function enhance = time_enhancing(t,thalf,h)

scale = 2^(-h); 
enhance = scale.*(1+tanh(3./thalf.*(t-thalf))).^(h);

end

%Function that influences the damepning of filopodia birth
function damp = time_dampening(t,p)
    damp = polyval(p,t);
    damp = max(1e-4,damp);
end

function f = feedback(B50,B)
    f = B50./(B50 + B);
end

function [c1_sF,c1_LF,c2_sF,c2_LF,c3, B50, c4, c5, p, thalf, h] = getParameters(param)

load('AllParameters');

c1_sF = Parameters.(char(param)).c1_sF;%birth           
c1_LF = Parameters.(char(param)).c1_LF;%birth
c2_sF = Parameters.(char(param)).c2_sF;%death
c2_LF = Parameters.(char(param)).c2_LF;%death
%bulbous dynamics
c3 = Parameters.(char(param)).c3;%rF2B
B50 = Parameters.(char(param)).B50;%rB2S
c4 = Parameters.(char(param)).c4;%rB20
c5 = Parameters.(char(param)).c5;%rB2S
% time dependent functions
% a) dampening of filo dynamics
p = Parameters.(char(param)).p;
% b) increase in propensity to form bulbous
thalf = Parameters.(char(param)).thalf;
h = Parameters.(char(param)).h;


end

  