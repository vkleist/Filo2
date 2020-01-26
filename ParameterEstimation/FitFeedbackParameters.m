close all
clear all

%Fits parameters c3, c5 and B50 from Table S3 

Mutants = {'Control';'Atg6_mutant'; 'Atg7_mutant';'Atg6_rescue'};
%choose mutant from the list above
mutant = 3;
MutantName = strrep(char(Mutants(mutant)),'_','-');
%danger: there is actually no real fit for trio
Colors = [.5 .5 .5;0 0 1;1 0 0;0 1 0];
Ref_yellow = [0.9 0.8 0.2];%

FSP = 20; % finite state projection
n = 6;

c = Colors(mutant,:);

Mutants(mutant)

%Expected numbers from PlotFiloNumbers.m
%Filopodia dynamics @P60
             %Control   Atg6_mutant     Atg7_mutant     Atg6_rescue'
sF          = [2.6      2               2.2             2.9];% average number short lived Filos @P60
LF          = [2.6      4.8             5.6             1.4];% average number long lived Filos @P60
d_sF        = [0.4348   0.3704          0.4167          0.4545];
d_LF        = [0.0667   0.0476          0.05            0.0769];

%
t = 20*60;%P60

% Parameters for slow time scale dynamics, previously determined in 
% Ozel (2019) Dev Cell, https://doi.org/10.1016/j.devcel.2019.06.014
halfmax = 1000;
exponent = 1;
enhance = enhancing(t,halfmax,exponent);
p = [-2.9740e-17   3.3115e-13  -1.2896e-09   2.0637e-06  -1.4483e-03   1.0021e+00];
damp = dampening(p,t);

%expected number of sF
Nbrs_sF = sF(mutant);%
%expected number of LF
Nbrs_LF = LF(mutant);%
F = Nbrs_sF+Nbrs_LF;

birth_sF = sF.*d_sF./damp
birth_LF = LF.*d_LF./damp

StartParamFB = log([0.022*enhance*F 0.1133 0.0282]);
%------- Parameters c4 & c6; see supplementary material of Kiral (2020) Nature Communications

%                 Control   Atg6_mutant     Atg7_mutant     Atg6_rescue 
DeathRate_sB    = [1/120    1/120           1/120           0.0712];%%death rate of short-lived bulbous; c4a
DeathRate_synB  = [1/133    1/133           1/133           1/133]; %death rate of long-lived bulbous; c4b


d_sB = DeathRate_sB(mutant);
d_synB = DeathRate_synB(mutant);

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
                      
%
edges = 0:7;
x_pos = 0:n;
Data_synBulb = ExpDensity_synBulb(mutant,:);
Data_sBulb = ExpDensity_sBulb(mutant,:);
Data = [Data_sBulb , Data_synBulb];

figure(1)
hold on
bar(edges(1:end-1),ExpDensity_synBulb(mutant,:)','EdgeColor','none','FaceColor',c);
hold on

figure(2)
hold on
bar(edges(1:end-1),ExpDensity_sBulb(mutant,:)','EdgeColor','none','FaceColor',c);
hold on

figure(3)
hold on
bar(edges(1:end-1),ExpDensity_AllBulb(mutant,:)','EdgeColor','none','FaceColor',c);
hold on

%% no feedback
StartParam = [rand rand];
options = optimset('Display','notify');%'iter');
x = fminsearch(@fitfun2,StartParam,options,Data,d_sB,d_synB,n);

c3 = exp(x(1));
c5 = c3/2*(1+sin(x(2)));% 0< c5 < c3;   %exp(x(2));
lambda1 = c3/(d_sB + c5);

%first distribution
p1 = poisspdf(x_pos,lambda1);

%exp(x(2));
lambda2 = c5./(d_synB).*lambda1;%sum((x_pos.*p1));%lambda1*(1+c5/d_synB);
p2 = poisspdf(x_pos,lambda2);

disp('Prob. Bulb becomes synaptogene (no feedback):')
PsB2synB = c5/(c5+d_sB)


figure(1)
stairs(x_pos-0.5,p2,'k:','Linewidth',4)
%bar(edges(1:end-1),[ExpDensity(mutant,:)',P]);
ylabel('probability','FontSize',14)
xlabel('number synaptogenic bulbous bulbous tips/(GC and time)','FontSize',14)
title(strcat(MutantName,': synaptogenic bulbs'),'FontSize',16)
%print(1,'-dpdf',strcat('./Figures/NoFeedback',char(Mutants(mutant)),'.pdf'))

figure(2)
stairs(x_pos-0.5,p1,'k:','Linewidth',4)
%bar(edges(1:end-1),[ExpDensity(mutant,:)',P]);
ylabel('probability','FontSize',14)
xlabel('number short-lived bulbous tips/(GC and time)','FontSize',14)
title(strcat(MutantName,': short-lived bulbs'),'FontSize',16)

Conv_PDF = fliplr(p1'*p2);
statdistr = nan(1,n+1);
for i = 0:n
    statdistr(i+1) = sum(diag(Conv_PDF,n-i));
end

figure(3)
ylabel('probability','FontSize',14)
xlabel('number bulbous tips/(GC and time)','FontSize',14)
title(strcat(MutantName,': all bulbs'),'FontSize',16)

%% Data fitting feedback (product inhibition on f1 or f2)
if mutant ==1
    
    StartParam = StartParamFB;%[r,rand,rand];
else
    r = rand;
    StartParam = StartParamFB;%[r,r-(1-sqrt(1+rand^2)),rand];
end
options = optimset('Display','iter');
x = fminsearch(@fitfun,StartParam,options,Data,FSP,d_sB,d_synB,n,mutant);


if mutant == 1 
    disp('Feedback on f1')
    c3_P60 = exp(x(1));%
 %   convert c3_P60 to c3
    c3 = c3_P60./(F.*enhance)
    c5 = exp(x(2))%
    B50 = exp(x(3))%

else 
    
    c3_P60 = exp(x(1));%
 %   convert c3_P60 to c3
    c3 = c3_P60./(F.*enhance)
    c5 = exp(x(2))%
    B50 = exp(x(3))%

end
 
L = makeGenerator(d_sB,d_synB,c3_P60,B50,FSP,c5);

disp('Prob. Bulb becomes synaptogene (feedback model, feedback off):')
PsB2synB = c5/(c5+d_sB)

maxIDx = length(Data);

M = L';
[V,D] = eig(M);
idx = find(round(diag(D).*1e5)./1e5 == 0);
statdistr = V(:,idx)./sum(V(:,idx));
TupelM = Idx2Tupel(FSP);
%add together sBulb + synBulb
P1 = zeros(n+1,1);
P2 = zeros(n+1,1);
for j = 0:n
    idx1 = TupelM(:,1) == j;
    P1(j+1) = sum(statdistr(idx1));
    idx2 = TupelM(:,2) == j;
    P2(j+1) = sum(statdistr(idx2));
end

%%
Expect_f = 0;
for j = 0:n
    Expect_f = Expect_f + P2(j+1) * feedback(B50,j);
end
disp('average r_5 @P60')
if mutant == 5
    (x_pos*P1)*c5.*Expect_f
    Expect_f = 1;
else
    (x_pos*P1)*c5
end

disp('average feedback f1 @P60')
Expect_f
disp('average r_2B @P60')
Expect_r2B = c3_P60
disp('average r_3 @P60')
c3_P60*Expect_f
disp('average r_4 @P60')
(x_pos*P1)*d_sB
disp('Average bulbs @P60 (measured)')
x_pos*ExpDensity_AllBulb(mutant,:)'

%

%prob sB
idx = TupelM(:,1)>0;
Prob_sB = sum(statdistr(idx));

%prob synB
idx = TupelM(:,2)>0;
Prob_synB = sum(statdistr(idx));

figure(1)
hold on
stairs(x_pos-0.5,P2,'-','Linewidth',4,'Color',[0.2 0.2 0.2])
legend('data','no feedback', 'feedback')
Path = './FiguresRobin/';
Name = strcat('Fig_BulbousFit_Synaptogenic',MutantName,'.tiff');
%print(1,'-dtiff',strcat(Path,Name,'.tiff'))

figure(2)
hold on
%stairs(x_pos-0.5,P1,'--','Linewidth',4,'Color','r')
stairs(x_pos-0.5,P1,'-','Linewidth',4,'Color',[0.2 0.2 0.2])
%if mutant == 1
    legend('data','no feedback', 'feedback')
%end
Path = './FiguresRobin/';
Name = strcat('Fig_BulbousFit_Transient',MutantName,'.tiff');
%print(2,'-dtiff',strcat(Path,Name,'.tiff'))

Prob_allB = nan(n+1,1);
%prob all Bulbs
for j = 0:n
    idx = sum(TupelM,2)==j;
    Prob_allB(j+1) = sum(statdistr(idx));
end
figure(3)
stairs(x_pos-0.5,Prob_allB,'-','Linewidth',5,'Color',[0.2 0.2 0.2])
%bar(edges(1:end-1),[ExpDensity(mutant,:)',P]);
set(gca,'FontSize',18)
ylabel('probability','FontSize',20)
xlabel('number bulbous tips/(GC and time)','FontSize',20)
title(strcat(MutantName,': all bulbs'),'FontSize',20)
%if mutant == 1
    legend('data','fit','FontSize',18);%, 'feedback')
% % %end
Path = './FiguresRobin/';
Name = strcat('Fig_BulbousFit_All',MutantName);
%print(3,'-dtiff',strcat(Path,Name,'.tiff'))
%print(3,'-depsc2',strcat(Path,Name,'.eps'))


disp('Average bulbs @P60 (predicted)')
x_pos*Prob_allB


function res = fitfun(x,Data,FSP,d_sB,d_synB,n,mutant)
    if mutant == 1 
        c3 = exp(x(1));%
        c5 = exp(x(2));%
        B50 = exp(x(3));%
    else 
        c3 = exp(x(1));%
        c5 = exp(x(2));%c3-(1-sqrt(1+x(2)^2));%c3/2*(1+sin(x(2)));%
        B50 = exp(x(3));%
    end
    L = makeGenerator(d_sB,d_synB,c3,B50,FSP,c5);

    M = L';
    [V,D] = eig(M);
    idx = find(round(diag(D).*1e5)./1e5 == 0);
    statdistr = V(:,idx)./sum(V(:,idx));
    
    % add together sBulb + synBulb
    TupelM = Idx2Tupel(FSP);
    % add together sBulb + synBulb
    P1 = zeros(1,n+1);
    P2 = zeros(1,n+1);
    for j = 0:n
        idx1 = TupelM(:,1) == j;
        P1(j+1) = sum(statdistr(idx1));
        idx2 = TupelM(:,2) == j;
        P2(j+1) = sum(statdistr(idx2));
    end
    
    res = KL(Data,[P1, P2]);
end

function L = makeGenerator(d_sB,d_synB,c3,B50,FSP,c5)
    L = zeros((FSP+1)^2,(FSP+1)^2);
    
    for i = 0:FSP % number short lived bulb
        for j = 0:FSP %number long lived buld
            idx = i*(FSP+1)+(j+1);%current state
            
            if i > 0 %sB -> *
                idx_death_i = (i-1)*(FSP+1)+(j+1);%(i-1,j)
                L(idx,idx_death_i) = d_sB*i;
            end
            
            if j > 0 %synB -> *
                idx_death_j = (i)*(FSP+1)+(j);%(i,j-1)
                L(idx,idx_death_j) = (d_synB)*j;
            end
            
            if i < FSP %F -> sB +1
                idx_birth_i = (i+1)*(FSP+1)+(j+1);%(i+1,j)
                L(idx,idx_birth_i) = c3.*feedback(B50,j);%c3.*feedback(B50,(i+j)); %c3.*feedback(B50,j);%
            end
            
            if j < FSP && i > 0% sB -> synB
                idx_birth_j = (i-1)*(FSP+1)+(j+2);%(i-1,j+1)
                L(idx,idx_birth_j) = (i)*c5;%
            end

        end
    end
    L = L - diag(sum(L,2));
end

function res = fitfun2(x,Data,d_sB,d_synB,n)
    %transformation [parameters have to be > 0]
    c3 = exp(x(1));
    c5 = c3/2*(1+sin(x(2)));
    
    %first distribution
    x_pos = 0:n;
    lambda1 = c3/(d_sB + c5);
    p1 = poisspdf(x_pos,lambda1);
    
    %exp(x(2));
    lambda2 = sum(c5.*p1./(d_synB));%lambda1*(1+c5/d_synB);
    
    p2 = poisspdf(x_pos,lambda2);
    statdistr = [p1,p2];
    res = KL(Data,statdistr);
end


function f = feedback(B50,B)
    f = B50./(B50 + B);
end

function p = poisspdf(x,lambda)
    p = lambda.^x./factorial(x).*exp(-lambda); 
end

function TupelM = Idx2Tupel(n)
    TupelM = nan((n+1)^2,2);
    %This Matrix contains at row i the number of short-lived bulbs (first
    %column) and long lived bulbs (second column)
    counter = 0;
    for i = 0:n
        for j = 0:n
            counter = counter +1;
            TupelM(counter,1) = i;
            TupelM(counter,2) = j;
        end
    end
end

function enhance = enhancing(t,thalf,h)
scale = 2^(-h); 
enhance = scale.*(1+tanh(3./thalf.*(t-thalf))).^(h);
end

function d = dampening(p,t)
    d = polyval(p,t);
    d = max(1e-4,d);
end
