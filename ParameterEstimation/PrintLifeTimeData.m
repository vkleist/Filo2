close all
clear all

% print the lifetimes of filopodia (Supplementary Table S1)

load('AllData.mat');


Mutants = {'Control';'Atg6_mutant'; 'Atg7_mutant'; 'Atg6_rescue'};
types = {'F';'sF';'All'};

TmpP60 = [];
AllTmpP60 = [];
edges1 = 0:1:60;
edges2 = 0:3:60;

Var1 = cell(length(Mutants)-1,1);
Var2 = cell(length(Mutants)-1,1);
Var3 = cell(length(Mutants)-1,1);
for i = 1:length(Mutants)
    TmpP60 = [];
    for j = 1:length(types)-1
        type = char(types(j));
        
        TmpP60 = [TmpP60;Output.(char(Mutants(i))).P60.(type).LTimes];
          
        m = mean(Output.(char(Mutants(i))).P60.(type).LTimes);
        s = std(Output.(char(Mutants(i))).P60.(type).LTimes);
        
        if j == 1
            Var1{i} = strcat(num2str(m,2),'(',num2str(s,2),')');
        else
            Var2{i} = strcat(num2str(m,2),'(',num2str(s,2),')');
        end
    end
    m = mean(TmpP60);
    s = std(TmpP60);
    Var3{i} = strcat(num2str(m,2),'(',num2str(s,2),')');
 
end


Names = {'Mutant';'F';'sF';'All'};
%print data to table
T = table(Mutants,Var1,Var2,Var3);%,'RowNames',Mutants);
T.Properties.VariableNames = Names;

% Now use this table as input in our input struct:
% LaTex table caption:
input.tableCaption = 'Average (standard deviation) life time of filopodia (min)';
input.data = T;
% Switch transposing/pivoting your table if needed:
input.transposeTable = 0;
% Switch to generate a complete LaTex document or just a table:
input.makeCompleteLatexDocument = 0;
% Now call the function to generate LaTex code:
latex = latexTable(input);    



