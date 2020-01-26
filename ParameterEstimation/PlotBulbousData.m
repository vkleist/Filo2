clear all
close all

load('AllData.mat');

Mutants = {'Control';'Atg6_mutant'; 'Atg7_mutant';'Atg6_rescue'};

LifeTimeThreshold = 40;%min

Var1 = cell(length(Mutants),1); % mean (std) nr bulbous tips

Var2 = cell(length(Mutants),1); % prob that there are 0,1,...5
Var3 = cell(length(Mutants),1); % prob that there are 0,1,...5
Var4 = cell(length(Mutants),1); % prob that there are 0,1,...5
Var5 = cell(length(Mutants),1); % prob that there are 0,1,...5
Var6 = cell(length(Mutants),1); % prob that there are 0,1,...5
Var7 = cell(length(Mutants),1); % prob that there are 0,1,...5
Var8 = cell(length(Mutants),1); % prob that there are 0,1,...5


Var2_1 = cell(length(Mutants),1); % mean (std) nr bulbous tips
Var2_2 = cell(length(Mutants),1); % prob that there are 0,1,...5
Var2_3 = cell(length(Mutants),1); % prob that there are 0,1,...5
Var2_4 = cell(length(Mutants),1); % prob that there are 0,1,...5
Var2_5 = cell(length(Mutants),1); % prob that there are 0,1,...5
Var2_6 = cell(length(Mutants),1); % prob that there are 0,1,...5
Var2_7 = cell(length(Mutants),1); % prob that there are 0,1,...5
Var2_8 = cell(length(Mutants),1); % prob that there are 0,1,...5

%Compute Nr BulbousTips per Time Instance

ReSample = 0;%100;%for each GrowthCone; if 0 then no resample
MaxNr = 6;
PSynapse = zeros(1,length(Mutants));

for j = 1:length(Mutants)
    mutant = char(Mutants(j));
    MutantName = strrep(mutant,'_','-');
    GCs = unique((Output.(mutant).P60.Bulb.GC));
    NrGCs = length(GCs);
    
    BulbTips = Output.(mutant).P60.Bulb.LTimes;
    %--------------
    lgx = BulbTips>=LifeTimeThreshold;
    
    Long = lgx;
    Short = ~lgx;

    Nbrs_sBulbous = zeros(NrGCs*1,60);
    Nbrs_synBulbous = zeros(NrGCs*1,60);
    ReSampleCounter = 0;
    for counter1 = GCs'
        %long lived
        Idx = Output.(mutant).P60.Bulb.GC == counter1 & Long;
        NrSynBulbous = sum(Idx);
        Tmp = Output.(mutant).P60.Bulb.StartTimes+1;
        StartsSynBF =Tmp(Idx); %Bulbous Filo
        Tmp = Output.(mutant).P60.Bulb.LTimes;
        LT_SynBT = Tmp(Idx);%LifeTime of Bulbous Tip
        Tmp = Output.(mutant).P60.Bulb.EndTimes+1;
        Ends_syn = Tmp(Idx);%Bulbous Tip
        
        %short-lived
        Idx = Output.(mutant).P60.Bulb.GC == counter1 & Short;
        NrSBulbous = sum(Idx);
        Tmp = Output.(mutant).P60.Bulb.StartTimes+1;
        StartsSBF =Tmp(Idx); %Bulbous Filo
        Tmp = Output.(mutant).P60.Bulb.LTimes;
        LT_SBT = Tmp(Idx);%LifeTime of Bulbous Tip
        Tmp = Output.(mutant).P60.Bulb.EndTimes+1;
        Ends_s = Tmp(Idx);%Bulbous Tip
        
        Starts_syn = StartsSynBF';
        Starts_s = StartsSBF';
         
            
        AllSynBulbous = zeros(1,60);
        for counter2 = 1:NrSynBulbous
            StartT = Starts_syn(counter2);
            EndT = Ends_syn(counter2);
            AllSynBulbous(StartT:EndT) = AllSynBulbous(StartT:EndT)+1;
        end

        AllSBulbous = zeros(1,60);
        for counter2 = 1:NrSBulbous
            StartT = Starts_s(counter2);
            EndT = Ends_s(counter2);
            AllSBulbous(StartT:EndT) = AllSBulbous(StartT:EndT)+1;
        end
        %plot Number of Filopodia at time instance
        figure(100+j)
        hold on
        plot(AllSynBulbous,':','LineWidth',2)

        figure(300+j)
        hold on
        plot(AllSBulbous,':','LineWidth',2)

        ReSampleCounter = ReSampleCounter+1;
        Nbrs_synBulbous(ReSampleCounter,:) = Nbrs_synBulbous(ReSampleCounter,:) + AllSynBulbous;
        Nbrs_sBulbous(ReSampleCounter,:) = Nbrs_sBulbous(ReSampleCounter,:) + AllSBulbous;

    end % end over growth cone
    %short lived filopodia
    data_syn = Nbrs_synBulbous(:);
    m_syn = mean(data_syn);
    s_syn = std(data_syn);
    
   %save data 
    Var1{j} = strcat(num2str(m_syn,2),'(',num2str(s_syn,2),')'); 
    
    data_s = Nbrs_sBulbous(:);
    m_s = mean(data_s);
    s_s = std(data_s);
    %save data 
    Var2_1{j} = strcat(num2str(m_syn,2),'(',num2str(s_syn,2),')'); 
    
    %Plot mean number of bulbous at time instance
    figure(100+j)
    hold on
    plot(mean(Nbrs_synBulbous),'k-','LineWidth',3)
    title(strcat('Nr. synaptogenic bulbous tips (',MutantName,')'),'FontSize',14)
    ylabel('P60','FontWeight','bold','FontSize',14)
    xlabel('time (min)','FontWeight','bold','FontSize',12)
    ylim([0 5])
    print(100+j,'-depsc2',strcat('./Figures/SynBulbousNumbers',MutantName,'.eps'))
    
    figure(300+j)
    hold on
    plot(mean(Nbrs_sBulbous),'k-','LineWidth',3)
    title(strcat('Nr. short-lived bulbous tips (',MutantName,')'),'FontSize',14)
    ylabel('P60','FontWeight','bold','FontSize',14)
    xlabel('time (min)','FontWeight','bold','FontSize',12)
    ylim([0 5])
    print(300+j,'-depsc2',strcat('./Figures/ShortBulbousNumbers',MutantName,'.eps'))
    
    %Plot number histogram
    edges = (0:7)-0.5;
    figure(200+j)
    hold on
    title(strcat('Nr. synaptogenic bulbous  tips (',MutantName,')'),'FontSize',14)
    [counts,edges] = histcounts(data_syn,edges);
    X = counts./sum(counts);
    histogram('BinEdges',edges,'BinCounts',X)
    line([m_syn  m_syn],[0 0.2],'Color','r','LineStyle','--','LineWidth',3)
    line([m_syn+s_syn m_syn+s_syn],[0 0.2],'Color','r','LineStyle',':','LineWidth',2)
    line([m_syn-s_syn  m_syn-s_syn],[0 0.2],'Color','r','LineStyle',':','LineWidth',2)
    ylim([0 1])
    ylabel('Frequency/Probability','FontWeight','bold','FontSize',14)
    xlabel('Number bulbous tips')
    print(200+j,'-depsc2',strcat('./Figures/SynBulbousNumberDistribution',MutantName,'.eps'));

    Var2{j} = strcat(num2str(X(1),2)); 
    Var3{j} = strcat(num2str(X(2),2)); 
    Var4{j} = strcat(num2str(X(3),2)); 
    Var5{j} = strcat(num2str(X(4),2)); 
    Var6{j} = strcat(num2str(X(5),2)); 
    Var7{j} = strcat(num2str(X(6),2)); 
    Var8{j} = strcat(num2str(X(7),2)); 
    
    
    %Plot number histogram
    edges = (0:7)-0.5;
    figure(400+j)
    hold on
    title(strcat('Nr. short-lived bulbous  tips (',MutantName,')'),'FontSize',14)
    [counts,edges] = histcounts(data_s,edges);
    X = counts./sum(counts);
    histogram('BinEdges',edges,'BinCounts',X)
    line([m_s  m_s],[0 0.2],'Color','r','LineStyle','--','LineWidth',3)
    line([m_s+s_s m_s+s_s],[0 0.2],'Color','r','LineStyle',':','LineWidth',2)
    line([m_s-s_s  m_s-s_s],[0 0.2],'Color','r','LineStyle',':','LineWidth',2)
    ylim([0 1])
    ylabel('Frequency/Probability','FontWeight','bold','FontSize',14)
    xlabel('Number bulbous tips')
    print(400+j,'-depsc2',strcat('./Figures/ShortBulbousNumberDistribution',MutantName,'.eps'));

    Var2_2{j} = strcat(num2str(X(1),2)); 
    Var2_3{j} = strcat(num2str(X(2),2)); 
    Var2_4{j} = strcat(num2str(X(3),2)); 
    Var2_5{j} = strcat(num2str(X(4),2)); 
    Var2_6{j} = strcat(num2str(X(5),2)); 
    Var2_7{j} = strcat(num2str(X(6),2)); 
    Var2_8{j} = strcat(num2str(X(7),2)); 
    
end
        
Names = {'Mutant';'Nr_syn_Bulbous'};
%print data to table
T = table(Mutants,Var1);%,'RowNames',Mutants);
T.Properties.VariableNames = Names;

% Now use this table as input in our input struct:
% LaTex table caption:
input.tableCaption = 'Number (standard deviation) of synaptogenic bulbous tips per growth cone (P60) and time instance';
input.data = T;
% Switch transposing/pivoting your table if needed:
input.transposeTable = 0;
% Switch to generate a complete LaTex document or just a table:
input.makeCompleteLatexDocument = 0;
% Now call the function to generate LaTex code:
latex = latexTable(input);    

Names = {'Mutant';'Nr_s_Bulbous'};
%print data to table
T = table(Mutants,Var2);%,'RowNames',Mutants);
T.Properties.VariableNames = Names;

% Now use this table as input in our input struct:
% LaTex table caption:
input.tableCaption = 'Number (standard deviation) of short-lived bulbous tips per growth cone (P60) and time instance';
input.data = T;
% Switch transposing/pivoting your table if needed:
input.transposeTable = 0;
% Switch to generate a complete LaTex document or just a table:
input.makeCompleteLatexDocument = 0;
% Now call the function to generate LaTex code:
latex = latexTable(input);    


%%

Names = {'Mutant';'n0';'n1';'n2';'n3';'n4';'n5';'n6'};
%print data to table
T = table(Mutants,Var2,Var3,Var4,Var5,Var6,Var7,Var8);%,'RowNames',Mutants);
T.Properties.VariableNames = Names;

% Now use this table as input in our input struct:
% LaTex table caption:
input.tableCaption = 'Probability the the number of synaptogenic bulbous tips per growth cone (P60) and time instance is n';
input.data = T;
% Switch transposing/pivoting your table if needed:
input.transposeTable = 0;
% Switch to generate a complete LaTex document or just a table:
input.makeCompleteLatexDocument = 0;
% Now call the function to generate LaTex code:
latex = latexTable(input);  

Names = {'Mutant';'n0';'n1';'n2';'n3';'n4';'n5';'n6'};
%print data to table
T = table(Mutants,Var2_2,Var2_3,Var2_4,Var2_5,Var2_6,Var2_7,Var2_8);%,'RowNames',Mutants);
T.Properties.VariableNames = Names;

% Now use this table as input in our input struct:
% LaTex table caption:
input.tableCaption = 'Probability the the number of short-lived bulbous tips per growth cone (P60) and time instance is n';
input.data = T;
% Switch transposing/pivoting your table if needed:
input.transposeTable = 0;
% Switch to generate a complete LaTex document or just a table:
input.makeCompleteLatexDocument = 0;
% Now call the function to generate LaTex code:
latex = latexTable(input);  


