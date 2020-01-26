%
%Implements and simulates (Gillespie algorithm) the mechanistic model (originally developed in
%Ozel et al (2019) Dev Cell) and described in the supplementary methods of
%Kiral et al (2020) Nature Communications that models a competition
%between filopodia for a limiting resource that stabilizes the filopodium.


close all
clear all

try
    load('AllSimulations.mat')
catch
end

Mutants = {'Control';'Atg6_mutant'; 'Atg7_mutant';'Atg6_rescue'};
mutantNr = 4;
mutant = char(Mutants(mutantNr));
MutantName = strrep(char(Mutants(mutantNr)),'_','-');
Colors = [.5 .5 .5;0 0 1;1 0 0;1 0 1];

if mutantNr == 4
    CompAdv = false;
    c0 = 0.0712; % from the parameter table
else
    CompAdv = true;%
    c0 = 1;%
end
ResourceLimitation = true;


% from Test_Feedback...
                %'WT';  'atg6'; 'atg7';     'atg6_rescue'
MutantFactorC3 = ones(1,4);%[1      0.98    2.09       1.145];%1.084]; % relative change in r2B
                %'WT';  'atg6'; 'atg7';     'atg6_rescue'
AverageBulbs =  [1.65    3.03    2.5         1.64];%0.45];%1.1210];

ResourceAccessFactor = AverageBulbs;

c = Colors(mutantNr,:);

%created by "PlotBulbousData.m"
                    % n0    n1      n2      n3       n4     n5      n6
ExpDensity_synBulb = [0     0.58    0.42    0        0      0       0; %Control
                      0     0.2     0.26    0.15     0.38   0       0; %Atg6_mutant
                      0     0.17    0.42    0.42     0      0       0; %Atg7_mutant  
                      0.76  0.24       0       0        0      0       0]; %Atg6_rescue 
        
                    % n0    n1      n2      n3       n4     n5      n6
ExpDensity_sBulb = [0.79   0.19     0.021   0       0       0       0; %Control
                    0.7    0.3      0       0       0       0       0; %Atg6_mutant
                    0.75   0.25     0       0       0       0       0; %Atg7_mutant  
                    0.13   0.47     0.27    0.11    0.012   0.0071  0]; %Atg6_rescue   
                    
                % if LifeTimeThreshold = 0 in counting routine
                    % n0    n1      n2      n3       n4     n5      n6
ExpDensity_AllBulb = [0     0.42    0.51    0.071    0       0       0; %Control                                               
                      0     0.083   0.28    0.25     0.3    0.087    0;  %Atg6_mutant                                          
                      0     0.047   0.43    0.49     0.031   0       0;  %Atg7_mutant
                      0.076 0.37    0.4     0.13     0.012   0.0071   0]; %Atg6_rescue
                                
% plots the data
                  
figure(3) 
edges = 0:7;
bar(edges(1:end-1),ExpDensity_AllBulb(mutantNr,:)','EdgeColor','none','FaceColor',c);
hold on

figure(4) 
edges = 0:7;
bar(edges(1:end-1),ExpDensity_synBulb(mutantNr,:)','EdgeColor','none','FaceColor',c);
hold on

figure(5) 
edges = 0:7;
bar(edges(1:end-1),ExpDensity_sBulb(mutantNr,:)','EdgeColor','none','FaceColor',c);
hold on

%parameters as described in the Suppl. Materials of Kiral et al. (2020) Nat
%Comm
r2B = 0.0948; % known (r2B wt)
c_in = 0.07; % 
c_out = 1.5;% 

NrStates = 120; % reflecting boundary

r2B = r2B*MutantFactorC3(mutantNr);

if ResourceLimitation
    X0 = zeros(1,NrStates+1);
    X0(1) = round(NrStates*ResourceAccessFactor(mutantNr));%R
end

BulbousLengthThreshold = round(NrStates/4); % below this bulbs are nor recogized as bulbs

t = 0;

X = X0;
tend = 100000;%one long trajectory

Xsave = nan(tend,length(X0));
tsave = nan(tend,1);
if ResourceLimitation
    FiloID = nan(X0(1)*2,3);  % lifetime, length, maximum length
else
    FiloID = nan(100,3);  % lifetime, length, maximum length
end
Lifetimes = nan(5000,4);  % Lifetime, Starttime, Endtime, max length
counter = 0;
counter2 = 0;
NextRecordingTime = t+1;%imaging frame interval (1 min)
while t < tend
    if CompAdv
        a = propensityCompetAdv(X,NrStates,r2B,c_in,c0,c_out,mutantNr);
    else
        a = propensityNoCompetAdv(X,NrStates,r2B,c_in,c0,c_out,mutantNr);
    end
    
    a0 = sum(a);
    
    r1 = rand;
    tau = (1/a0)*log(1/r1);%time step
    t = t+tau; %time update
    % 2. sample reaction
    r2 = rand;
    j = find(r2 <= cumsum(a)/a0,1);
    
    %save data
   if t >= NextRecordingTime
        CurrentRecordingTime = floor(t);
        NrRecordings = CurrentRecordingTime-NextRecordingTime +1;
        
        Xsave(NextRecordingTime+1:CurrentRecordingTime+1,:) = repmat(X,NrRecordings,1);
        tsave(NextRecordingTime+1:CurrentRecordingTime+1) = repmat(t,NrRecordings,1);    
        
        NextRecordingTime = CurrentRecordingTime+1;
    end
    
    %ExecuteReaction
    [X,Lifetimes,FiloID,counter2] = ExecuteReaction(X,j,NrStates,Lifetimes,FiloID,tau,counter2,t);
    if ~ResourceLimitation
        X(1) = X0(1);
    end
end

Lifetimes = Lifetimes(1:counter2,:);

LT = Lifetimes(:,1);
idx = Lifetimes(:,4) > BulbousLengthThreshold;

MechSimulation.(mutant).P60.Bulb.LTimes = LT(idx);
save('./AllSimulations.mat','MechSimulation');%,'-append');


figure(6)
histogram(LT(idx),[0:10:50, inf],'Normalization','probability','FaceColor',Colors(mutantNr,:));
xlim([0 60])
title('histogram of lifetimes','Fontsize',20)
set(gca,'FontSize',18)
set(gca,'XTickLabels', {'0';'10';'20';'30';'40';'50';'>60'},'Fontsize',18)
xlabel('Lifetime (min)','Fontsize',20)
ylabel('Frequency','Fontsize',20)
ylim([0 0.9])
%print(6,'-dtiff',strcat('./Figures/Mechanistic_LifeTimes',MutantName,'.tiff'))
%print(6,'-depsc2',strcat('./Figures/Mechanistic_LifeTimes',MutantName,'.eps'))

%% Number short vs. long lived Bulbs
NbrsStabBulb = zeros(1,tend);
NbrsBulbs = zeros(1,tend);
        
%For each growth cone
%Number of Bulbpodia at time instance
Idx =  Lifetimes(:,4)> BulbousLengthThreshold & Lifetimes(:,1) < 40;
NrBulbs = sum(Idx);
Starts = ceil(Lifetimes(Idx,2));
Ends = floor(Lifetimes(Idx,3));
AllNrBulbs = zeros(1,tend);
for counter2 = 1:NrBulbs
    StartT = Starts(counter2);
    EndT = Ends(counter2);
    AllNrBulbs(StartT:EndT) = AllNrBulbs(StartT:EndT)+1;
end
% %Number of Bulbpodia at time instance
Idx = Lifetimes(:,4)> BulbousLengthThreshold & Lifetimes(:,1) >= 40;
NrStabBulbs = sum(Idx);
Starts = ceil(Lifetimes(Idx,2));
Ends = floor(Lifetimes(Idx,3));
AllStabBulbs = zeros(1,tend);   
for counter2 = 1:NrStabBulbs
    StartT = Starts(counter2);
    EndT = Ends(counter2);
    AllStabBulbs(StartT:EndT)= AllStabBulbs(StartT:EndT)+ 1;
end

%plot Number of Bulbpodia at time instance
figure(101)
subplot(2,1,1)
hold on
plot(AllNrBulbs,':','LineWidth',2)
ylabel('short lived')
subplot(2,1,2)
hold on
plot(AllStabBulbs,':','LineWidth',2)
ylabel('long lived')

figure(3)%All Bulbs
edges = 0:8;
[N,edges] = histcounts(AllStabBulbs+AllNrBulbs,edges);
stairs(edges(1:end-1)-0.5,N./sum(N),'k','LineWidth',5)
set(gca,'FontSize',18)
title(strcat('Number bulbs:',MutantName),'Fontsize',20)
xlabel('number','Fontsize',20)
ylabel('frequency','Fontsize',20)
SimDensity_AllBulb = N./sum(N)

%print(3,'-dtiff',strcat('./Figures/Mechanistic_Bulb',MutantName,'.tiff'))
%print(3,'-depsc2',strcat('./Figures/Mechanistic_Bulb',MutantName,'.eps'))

figure(4)%long-lived
[N,edges] = histcounts(AllStabBulbs,edges);
stairs(edges(1:end-1)-0.5,N./sum(N),'k','LineWidth',3)
title('Number stable bulbs')
xlabel('number','Fontsize',16)
ylabel('frequency','Fontsize',16)
set(gca,'FontSize',14)
SimDensity_StabBulb = N./sum(N)

figure(5)% short-lived
[N,edges] = histcounts(AllNrBulbs,edges);
stairs(edges(1:end-1)-0.5,N./sum(N),'k','LineWidth',3)
title('Number transient bulbs')
xlabel('number','Fontsize',16)
ylabel('frequency','Fontsize',16)
set(gca,'FontSize',14)
SimDensity_TransBulb = N./sum(N)

function a = propensityCompetAdv(X,NrStates,r2B,c_in,c0,c_out,mutantNr)
    
    a = zeros(NrStates*3+1,1);
    %1)%0 -> B; filopodial birth
    a(1) = r2B; 
    %2)%B(i)+ R -> B(i+1); accumulation of resource
    a(2:NrStates) = X(1).*X(2:end-1)*c_in;
    %B(N)+R -> B(N+1)
    a(NrStates+1) = 0; % boundary condition
    %3)B(i) -> 0 + i*R ; retraction and release of resource
    current = NrStates+1;
    a(current+1:current+1+NrStates-1) = X(2:end).*c0.*1./((1:NrStates));
    %4) B(i) -> B(i-1) + R; release of a single resource
    current = current+1+NrStates-1;
    % B(0) -> B(-1)  
    a(current+1) = 0; % boundary condition
    % B(i) -> B(i-1) + R
    a(current+2:current+1+NrStates-1) = X(3:end).*c_out; 

end

function a = propensityNoCompetAdv(X,NrStates,r2B,c_in,c0,c_out,mutantNr)
    r2B = r2B*1;
    a = zeros(NrStates*3+1,1);
    %1)0 -> B
    a(1) = r2B; 
    %2)B(i)+ R -> B(i+1)
    a(2:NrStates) = X(1).*X(2:end-1)*c_in;%
    a(NrStates+1) = 0;%B(N)+R -> B(N+1); boundary condition
    %3)B(i) -> 0 + i*R
    current = NrStates+1;
    a(current+1:current+1+NrStates-1) = X(2:end).*c0;% 
    %4) B(i) -> B(i-1) + R
    current = current+1+NrStates-1;
    a(current+1) = 0;% B(0) -> B(-1); boundary condition
    a(current+2:current+1+NrStates-1) = X(3:end).*c_out; % B(i) -> B(i-1) + R

end


function [X,Lifetimes,FiloID,counter2] = ExecuteReaction(X,r,NrStates,Lifetimes,FiloID,tau,counter2,t)
    %r is the index of the reaction
    
    idx = find(~isnan(FiloID(:,1)));
    FiloID(idx,1) = FiloID(idx,1) + tau;
    
    %1) 0 -> B
   if r == 1 
       X(2) = X(2) + 1;   
       idx = find(isnan(FiloID(:,1)),1);
       FiloID(idx,1) = 0;
       FiloID(idx,2) = 1;
       FiloID(idx,3) = 1;
   %2) R + B(i) -> B(i+1)
   elseif 2 <= r && r <= NrStates+1 
      X(r) = X(r)-1;
      X(r+1) = X(r+1) + 1;
      X(1) = X(1) - 1;
      
      leng = r-1;
      lgx = find(round(FiloID(:,2)) == leng);
      idx = randi(length(lgx));
      idx2 = lgx(idx);
      FiloID(idx2,2) = leng+1;
      FiloID(idx2,3) = max(FiloID(idx2,3),FiloID(idx2,2));
      
   %3) B(i) -> 0 + i*R
   elseif NrStates+2 <= r && r <= 2*NrStates+1 
       idx = r-NrStates;
       X(idx) = X(idx)-1;
       X(1) = X(1) + idx-2;
       
      leng = idx-1;
      lgx = find(round(FiloID(:,2)) == leng);
      idx = randi(length(lgx));
      idx2 = lgx(idx);
      
      if FiloID(idx2,1) >= t-floor(t)  % if lifetime is exceeds the last storage/imaging point
          counter2 = counter2 +1;
          Lifetimes(counter2,1) = FiloID(idx2,1);%LifeTime
          Lifetimes(counter2,3) = t;%EndTime
          Lifetimes(counter2,2) = t-Lifetimes(counter2,1);%StartTime
          Lifetimes(counter2,4) = FiloID(idx2,3);%max length
      end
      
      FiloID(idx2,3) = nan;
      FiloID(idx2,2) = nan;
      FiloID(idx2,1) = nan;

   %4) B(i) -> B(i-1) + R    
   else 
       idx = r-2*NrStates;
       X(idx-1) = X(idx-1)+1;
       X(idx) = X(idx) -1;
       X(1) = X(1) + 1;
       
       
      leng = idx-1;
      lgx = find(round(FiloID(:,2)) == leng);
      idx = randi(length(lgx));
      idx2 = lgx(idx);
      FiloID(idx2,2) = leng-1;
   end

end

