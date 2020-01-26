close all
clear all

%Loads the data and plots the number of filopodia with a superimposed
%Poisson distribution (Fig. S9)
%Prints the number and standard deviation of filopodia (Table S1)

load('AllData.mat');
Mutants = {'Control';'Atg6_mutant'; 'Atg7_mutant';'Atg6_rescue'};
Times = {'P60'};
Colors = [.5 .5 .5;0 0 1;1 0 0;0 1 0];

poisspdf = @(x,lambda) lambda.^x./factorial(x).*exp(-lambda); 
Var1 = cell(length(Mutants),1);
Var3 = cell(length(Mutants),1);

for j = 1:length(Mutants)
    mutant = char(Mutants(j));
    MutantName = strrep(mutant,'_','-');
    for i = 1:length(Times)    
        Exp = unique(Output.(mutant).(char(Times(i))).F.GC);
        NrExp = length(Exp);
        
        NbrsStabFilo = zeros(NrExp,60);
        NbrsFilos = zeros(NrExp,60);
        
        for counter1 = Exp'
            %For each growth cone
            %Number of Filopodia at time instance
            Idx = Output.(mutant).(char(Times(i))).F.GC == counter1;
            NrFilos = sum(Idx);
            Starts = Output.(mutant).(char(Times(i))).F.StartTimes(Idx)+1;
            Ends = Output.(mutant).(char(Times(i))).F.EndTimes(Idx)+1;
            AllNrFilos = zeros(1,60);
            for counter2 = 1:NrFilos
                StartT = Starts(counter2);
                EndT = Ends(counter2);
                AllNrFilos(StartT:EndT) = AllNrFilos(StartT:EndT)+1;
            end
            % %Number of stable Filopodia at time instance
            Idx = Output.(mutant).(char(Times(i))).sF.GC == counter1;
            NrStabFilos = sum(Idx);
            Starts = Output.(mutant).(char(Times(i))).sF.StartTimes(Idx)+1;
            Ends = Output.(mutant).(char(Times(i))).sF.EndTimes(Idx)+1;
            AllStabFilos = zeros(1,60);   
            for counter2 = 1:NrStabFilos
                StartT = Starts(counter2);
                EndT = Ends(counter2);
                AllStabFilos(StartT:EndT)= AllStabFilos(StartT:EndT)+ 1;
            end
            NbrsFilos(counter1,:) = NbrsFilos(counter1,:) + AllNrFilos;
            NbrsStabFilo(counter1,:) = NbrsStabFilo(counter1,:) + AllStabFilos;
        end % end Exp
        
        %short lived filopodia
        data = NbrsFilos(:);
        m = mean(data);
        s = std(data);
        %save data 
        if i == 1 % P60
            Var1{j} = strcat(num2str(m,2),'(',num2str(s,2),')');
        end  
              
        %Plot number histogram
        edges = (0:20)-0.5;
        figure(201+j)
        subplot(1,2,(i-1)*3+1)
        hold on
        title(strcat('sF (',MutantName,')'),'FontSize',18,'Interpreter','Latex')
        [counts,edges] = histcounts(data,edges);
        histogram('BinEdges',edges,'BinCounts',counts./sum(counts),'FaceColor',Colors(j,:))
        ylim([0 0.5])
        PoissPred = poisspdf(0:20,m);
        stairs(edges,PoissPred,'k-','LineWidth',5)
        line([m m],[0 0.35],'Color','k','LineWidth',4,'LineStyle',':')
        
        %long lived filopodia
        data = NbrsStabFilo(:);
        m = mean(data);
        s = std(data);

        %save data 
        if i == 1 % P60
            Var3{j} = strcat(num2str(m,2),'(',num2str(s,2),')');
        end  
        
        subplot(1,2,(i-1)*3+2)
        hold on
        title(strcat('$\ell$F (',MutantName,')'),'FontSize',18,'Interpreter','Latex')
        [counts,edges] = histcounts(data,edges);
        histogram('BinEdges',edges,'BinCounts',counts./sum(counts),'FaceColor',Colors(j,:))
        ylim([0 0.5])
        PoissPred = poisspdf(0:20,m);
        stairs(edges,PoissPred,'k-','LineWidth',5)
        line([m m],[0 0.35],'Color','k','LineWidth',4,'LineStyle',':')

    end % end over times i
    
    
    figure(201+j)
    subplot(1,2,1)
    set(gca,'FontSize',18)
    ylabel('Probability','FontWeight','bold','FontSize',20)
    xlabel('Number','FontWeight','bold','FontSize',20)
    subplot(1,2,2)
    set(gca,'FontSize',18)
    xlabel('Number','FontWeight','bold','FontSize',20)
    set(gcf,'renderer','Painters')
%    print(201+j,'-depsc2',strcat('./Figures/NumberDistributionVsPoisson',mutant,'.eps'));
%    print(201+j,'-dtiff',strcat('./Figures/NumberDistributionVsPoisson',mutant,'.tiff'));
end % end over Mutants

%Print mean numbers into table
Names = {'Mutant';'short_lived';'long_lived'};
%print data to table
T = table(Mutants,Var1,Var3);%,'RowNames',Mutants);
T.Properties.VariableNames = Names;

% Now use this table as input in our input struct:
% LaTex table caption:
input.tableCaption = 'Average (standard deviation) numbers of filopodia per time instance';
input.data = T;
% Switch transposing/pivoting your table if needed:
input.transposeTable = 0;
% Switch to generate a complete LaTex document or just a table:
input.makeCompleteLatexDocument = 0;
% Now call the function to generate LaTex code:
latex = latexTable(input);    