%Generates Figure 6h

close all
clear all


Mutants = {'Control';'Atg6_mutant'; 'Atg7_mutant';'Atg6_rescue'};
Colors = [.5 .5 .5;0 0 1;1 0 0;1 0 1];

XLabels = {};
counter = 0; 

for mutant = 1:length(Mutants)

    c = Colors(mutant,:);

    %created by "PlotBulbousData.m"
                        % n0    n1      n2      n3       n4     n5          n6
    ExpDensity_synBulb = [0     0.58    0.42    0        0      0           0; %Control
                          0     0.2     0.26    0.15     0.38   0           0; %Atg6_mutant
                          0     0.17    0.42    0.42     0      0           0; %Atg7_mutant  
                          0.76  0.24    0       0        0      0           0]; %Atg6_rescue 
    %created by "MechanisticModel.m"
    SimDensity_synBulb = [0.0678 0.4689 0.4431  0.0201   0.0001   0         0;
                          0.0116 0.0771 0.3194  0.5219   0.0673   0.0025    0.0001;  %Atg6_mutant    
                          0.0234 0.1529 0.5118  0.3021   0.0092   0.0006    0;  %Atg7_mutant     
                          0.7405 0.2249 0.0294  0.0045   0.0007   0.0000    0]; %Atg6_rescue
                        % n0    n1      n2      n3       n4     n5      n6
    ExpDensity_sBulb = [0.79   0.19     0.021   0       0       0       0; %Control
                        0.7    0.3      0       0       0       0       0; %Atg6_mutant
                        0.75   0.25     0       0       0       0       0; %Atg7_mutant  
                        0.13   0.47     0.27    0.11    0.012   0.0071  0]; %Atg6_rescue   
     
    SimDensity_sBulb = [0.8349 0.1541  0.0106  0.0004  0       0       0;%Control
                        0.8308 0.1500  0.0187  0.0005  0       0       0;  %Atg6_mutant    
                        0.8363 0.1525  0.0109  0.0003  0       0       0;  %Atg7_mutant     
                        0.3743 0.4037  0.1738  0.0417  0.0058  0.0007  0]; %Atg6_rescue

                        % n0    n1      n2      n3       n4     n5      n6
    ExpDensity_AllBulb = [0     0.42    0.51    0.071    0       0       0; %Control                                               
                          0     0.083   0.28    0.25     0.3    0.087    0;  %Atg6_mutant                                          
                          0     0.047   0.43    0.49     0.031   0       0;  %Atg7_mutant
                          0.076 0.37    0.4     0.13     0.012   0.0071   0]; %Atg6_rescue

    SimDensity_AllBulb = [0.0479 0.3584 0.5480  0.0450   0.0007  0        0;%Control
                          0.0084 0.0562 0.2478  0.5582   0.1203  0.0088   0.0003;  %Atg6_mutant    
                          0.0168 0.1155 0.4467  0.3967   0.0232  0.0012   0.0001;  %Atg7_mutant     
                          0.2633 0.3834 0.2594  0.0769   0.0147  0.0022   0.0001]; %Atg6_rescue
                 
   Mean_Exp = ExpDensity_AllBulb(mutant,:)*(0:length(ExpDensity_AllBulb(mutant,:))-1)';
   Stdev_Exp = sqrt(ExpDensity_AllBulb(mutant,:)*((0:length(ExpDensity_AllBulb(mutant,:))-1)'-Mean_Exp).^2);
   
   %All Bulbs
   figure(1)
   hold on
   bar((mutant-1)*2+1, Mean_Exp,'FaceColor',c)
   errorbar((mutant-1)*2+1,Mean_Exp,Stdev_Exp,'Color','k','LineWidth',1.5,'CapSize',15)
   counter = counter +1;
   XLabels{counter} = 'Experiment';
   
   Mean_Sim = SimDensity_AllBulb(mutant,:)*(0:length(SimDensity_AllBulb(mutant,:))-1)';
   Stdev_Sim = sqrt(SimDensity_AllBulb(mutant,:)*((0:length(SimDensity_AllBulb(mutant,:))-1)'-Mean_Sim).^2);
   bar((mutant-1)*2+2, Mean_Sim,'FaceColor',c)
   errorbar((mutant-1)*2+2,Mean_Sim,Stdev_Sim,'Color','k','LineWidth',1.5,'CapSize',15)
   counter = counter +1;
   XLabels{counter} = 'Model';
    
   % Synaptogenic
   Mean_Exp = ExpDensity_synBulb(mutant,:)*(0:length(ExpDensity_synBulb(mutant,:))-1)';
   Stdev_Exp = sqrt(ExpDensity_synBulb(mutant,:)*((0:length(ExpDensity_synBulb(mutant,:))-1)'-Mean_Exp).^2);
   figure(2)
   hold on
   bar((mutant-1)*2+1, Mean_Exp,'FaceColor',c)
   errorbar((mutant-1)*2+1,Mean_Exp,Stdev_Exp,'Color','k','LineWidth',1.5,'CapSize',15)
   counter = counter +1;
   XLabels{counter} = 'Experiment';
   
   Mean_Sim = SimDensity_synBulb(mutant,:)*(0:length(SimDensity_synBulb(mutant,:))-1)';
   Stdev_Sim = sqrt(SimDensity_synBulb(mutant,:)*((0:length(SimDensity_synBulb(mutant,:))-1)'-Mean_Sim).^2);
   bar((mutant-1)*2+2, Mean_Sim,'FaceColor',c)
   errorbar((mutant-1)*2+2,Mean_Sim,Stdev_Sim,'Color','k','LineWidth',1.5,'CapSize',15)
   counter = counter +1;
   XLabels{counter} = 'Model';
   
    % transient
   Mean_Exp = ExpDensity_sBulb(mutant,:)*(0:length(ExpDensity_sBulb(mutant,:))-1)';
   Stdev_Exp = sqrt(ExpDensity_sBulb(mutant,:)*((0:length(ExpDensity_sBulb(mutant,:))-1)'-Mean_Exp).^2);
   figure(3)
   hold on
   bar((mutant-1)*2+1, Mean_Exp,'FaceColor',c)
   errorbar((mutant-1)*2+1,Mean_Exp,Stdev_Exp,'Color','k','LineWidth',1.5,'CapSize',15)
   counter = counter +1;
   XLabels{counter} = 'Experiment';
   
   Mean_Sim = SimDensity_sBulb(mutant,:)*(0:length(SimDensity_sBulb(mutant,:))-1)';
   Stdev_Sim = sqrt(SimDensity_sBulb(mutant,:)*((0:length(SimDensity_sBulb(mutant,:))-1)'-Mean_Sim).^2);
   bar((mutant-1)*2+2, Mean_Sim,'FaceColor',c)
   errorbar((mutant-1)*2+2,Mean_Sim,Stdev_Sim,'Color','k','LineWidth',1.5,'CapSize',15)
   counter = counter +1;
   XLabels{counter} = 'Model';
   
   
end

figure(1)
set(gca,'FontSize',12);
xlim([0.5 8.5])
set(gca,'XTick',1:8,'XTickLabel',XLabels)
ylabel('Number/terminal')
title('All Bulbs')
set(get(gca,'xlabel'),'FontWeight','bold');
set(get(gca,'ylabel'),'Fontsize',16);
%print(1,'-dtiff',strcat('./Figures/Fig4_lowerPanelMechModel_AllBulbs.tiff'))
%print(1,'-depsc2',strcat('./Figures/Fig4_lowerPanelMechModel_AllBulbs.eps'))

figure(2)
set(gca,'FontSize',12);
xlim([0.5 8.5])
set(gca,'XTick',1:8,'XTickLabel',XLabels)
ylabel('Number/terminal')
title('synaptogenic Bulbs')
set(get(gca,'xlabel'),'FontWeight','bold');
set(get(gca,'ylabel'),'Fontsize',16);
%print(2,'-dtiff',strcat('./Figures/Fig4_lowerPanelMechModel_synaptoBulbs.tiff'))
%print(2,'-depsc2',strcat('./Figures/Fig4_lowerPanelMechModel_synaptoBulbs.eps'))
   

figure(3)
set(gca,'FontSize',12);
xlim([0.5 8.5])
set(gca,'XTick',1:8,'XTickLabel',XLabels)
ylabel('Number/terminal')
title('transient Bulbs')
set(get(gca,'xlabel'),'FontWeight','bold');
set(get(gca,'ylabel'),'Fontsize',16);
%print(3,'-dtiff',strcat('./Figures/Fig4_lowerPanelMechModel_transientBulbs.tiff'))
%print(3,'-depsc2',strcat('./Figures/Fig4_lowerPanelMechModel_transientBulbs.eps'))
   
