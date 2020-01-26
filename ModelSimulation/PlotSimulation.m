close all
clear all


Mutants = {'Control';'Atg6_mutant'; 'Atg7_mutant'; 'Atg6_rescue'};
MutantNames = {'Control';'Atg6-mutant'; 'Atg7-mutant'; 'Atg6-rescue'};

Colors = [.5 .5 .5;0 0 1;1 0 0;1 0 1];

XPositions = (1:100:3600);
XPositionsPlot = XPositions./60+40;

for i = 1:length(Mutants)
    mutant = char(Mutants(i));
    MutantName = char(MutantNames(i))
    Filename = strcat('./EnsembleData/DataSimple_',mutant,'.mat');
    S = load(Filename);
    PlotColor = Colors(i,:);
    
    %Filo
    Ref1 = S.ens_Data(:,XPositions,1);
    Ref2 = S.ens_Data(:,XPositions,2);
    MRef = mean(Ref1+Ref2);
    StdRef = std(Ref1+Ref2);
    figure(10)
    hold on
    plot(XPositionsPlot,MRef,'Color',PlotColor,'LineWidth',3)
    %plot(XPositionsPlot,MRef+StdRef,'Color',PlotColor,'LineWidth',1,'LineStyle',':')
    %plot(XPositionsPlot,MRef-StdRef,'Color',PlotColor,'LineWidth',1,'LineStyle',':')
    
    %Bulbs
    Ref1 = S.ens_Data(:,XPositions,3);
    Ref2 = S.ens_Data(:,XPositions,4);
    MRef = mean(Ref1+Ref2);
    StdRef = std(Ref1+Ref2);
    figure(8)
    hold on
    plot(XPositionsPlot,MRef,'Color',PlotColor,'LineWidth',3)
    %plot(XPositionsPlot,MRef+StdRef,'Color',PlotColor,'LineWidth',1,'LineStyle',':')
    %plot(XPositionsPlot,MRef-StdRef,'Color',PlotColor,'LineWidth',1,'LineStyle',':')

    %Synapses
    Ref1 = S.ens_Data(:,XPositions,5);
    MRef = mean(Ref1);
    StdRef = std(Ref1);
    figure(9)
    hold on
    plot(XPositionsPlot,MRef,'Color',PlotColor,'LineWidth',3)
    %plot(XPositionsPlot,MRef+StdRef,'Color',PlotColor,'LineWidth',1,'LineStyle',':')
    %plot(XPositionsPlot,MRef-StdRef,'Color',PlotColor,'LineWidth',1,'LineStyle',':')
%---------

    %All Filopodia


end

figure(10)
%xlabel('time (hours)')
%title(strcat('Number of filopodia (',AutophagyTypeName,')'))
title('Filopodia number (simulated)')
ylabel('Number/terminal')
set(get(gca,'xlabel'),'Fontsize',16);
set(get(gca,'ylabel'),'Fontsize',16);
set(get(gca,'title'),'Fontsize',16);
set(gca,'XTicklabel',cellstr(strcat('P',num2str((40:10:100)'))))
set(gca,'FontSize',14);
ylim([0 20])
%print(10,'-dtiff',strcat('./Figures/RidvanSimulationFilopodia.tiff'))


figure(8)
%xlabel('time (hours)')
title('Bulbous tip number (simulated)')
%title(strcat('Number of bulbous tips (',AutophagyTypeName,')'))
ylabel('Number/terminal')
set(get(gca,'xlabel'),'Fontsize',16);
set(get(gca,'ylabel'),'Fontsize',16);
set(get(gca,'title'),'Fontsize',16);
set(gca,'XTicklabel',cellstr(strcat('P',num2str((40:10:100)'))))
set(gca,'FontSize',14);
%set(gca,'FontWeight','b');
ylim([0 5])
%print(8,'-dtiff',strcat('./Figures/RidvanSimulationBulb.tiff'))


figure(9)
%xlabel('time (hours)')
%title(strcat('Number of synapses (',AutophagyTypeName,')'))
title('Synapses number (simulated)')
ylabel('Number/terminal')
set(get(gca,'xlabel'),'Fontsize',16);
set(get(gca,'ylabel'),'Fontsize',16);
set(get(gca,'title'),'Fontsize',16);
set(gca,'XTicklabel',cellstr(strcat('P',num2str((40:10:100)'))))
set(gca,'FontSize',14);
%set(gca,'FontWeight','b');
ylim([0 50])
%print(9,'-dtiff',strcat('./Figures/RidvanSimulationSynapse.tiff'))


