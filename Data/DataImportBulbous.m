close all
clear all


path = './Files/';
Mutants = {'Control';'Atg6_mutant'; 'Atg7_mutant';'Atg6_rescue'};
load('AllData.mat');


NrMutants = length(Mutants);
for j = 1:length(Mutants)
    %Path name and entering the path
    FullPath = strcat(path,strcat('Bulbous_',char(Mutants{j})));
    cd(FullPath)
    %check content and find files
    k = ls;
    FileTypePos = strfind(k,'.xlsx');
    NrFiles = length(FileTypePos);
    Filenames = cell(NrFiles,1);
    pos = 1;
    for i = 1:NrFiles
        filename = k(pos:FileTypePos(i)+4);
        Filenames{i} = filename;
        pos = FileTypePos(i)+6;
    end

        % save statistics for each mutant
        NrP60 = NrFiles;%sum(Assignment == 2); % Nr of P60 files
        BulbousInfoP60 = nan(20,5,NrP60);%For each Filo (dim1): Start, End, Lifetime (dim2), dim3 = growth cone, censored (1= lifetime could be larger)

        LifeTimeValues = nan(20,NrP60,1);
        NbrsBulbousAllExperimentsP60 = zeros(NrP60,60);
        counterP60 = 0;
        for z = 1:length(Filenames)
            filename = char(Filenames{z});
            GCPos = strfind(filename,'CG')+2;
            GC = str2num(filename(GCPos));
            %read table
            T = readtable(filename);
            
            NrEntries = length(T.BulbousTipID);

            StartTimes = T.BetweenFrames_x_y_;
            EndTimes = T.Var4;
            LifeTimes = T.LifeTime_min_;
            %save the data

            BulbousInfoP60(1:NrEntries,1,z) = StartTimes;
            BulbousInfoP60(1:NrEntries,2,z) = EndTimes;
            BulbousInfoP60(1:NrEntries,3,z) = LifeTimes;
            BulbousInfoP60(1:NrEntries,4,z) = GC;
            BulbousInfoP60(1:NrEntries,5,z) = (EndTimes == 59) | (StartTimes == 0); % right censored
            
            LifeTimeValues(1:length(LifeTimes),z) = LifeTimes;

             %% Compute a number of statistics
            %Plot Filo & Stabous Numbers
            BulbInfoTmp1lbs = BulbousInfoP60(:,:,z);

            AllNrBulbs = zeros(1,60);
            for counter = 1:NrEntries
                Start = BulbInfoTmp1lbs(counter,1)+1;
                End = BulbInfoTmp1lbs(counter,2)+1;
                AllNrBulbs(Start:End) = AllNrBulbs(Start:End)+1;
            end
            
            figure(101)
            subplot(1,NrMutants,j)
            hold on
            plot(AllNrBulbs,':','LineWidth',2)
            %------------------
            counterP60 = counterP60 +1;
            NbrsBulbousAllExperimentsP60(counterP60,:) = AllNrBulbs;
% 
        end %end over files
        % --- save data in data structure

        BulbousInfoP60_tmp = [];
        for i = 1:NrFiles
            BulbousInfoP60_tmp = [BulbousInfoP60_tmp;BulbousInfoP60(:,:,i)];
        end
        LGX = ~isnan(BulbousInfoP60_tmp(:,1));
        BulbousInfoP60_tmp2 = BulbousInfoP60_tmp(LGX,:);
        Output.(char(Mutants(j))).P60.Bulb.LTimes = BulbousInfoP60_tmp2(:,3);    
        Output.(char(Mutants(j))).P60.Bulb.StartTimes = BulbousInfoP60_tmp2(:,1); 
        Output.(char(Mutants(j))).P60.Bulb.EndTimes = BulbousInfoP60_tmp2(:,2); 
        Output.(char(Mutants(j))).P60.Bulb.GC = BulbousInfoP60_tmp2(:,4); 
        Output.(char(Mutants(j))).P60.Bulb.Censored = BulbousInfoP60_tmp2(:,5); 
        %-----
    cd ..
    cd ..
    
    figure(201)
    subplot(1,NrMutants,j)
    hold on
    edges = (0:7)-0.5;

    Data = NbrsBulbousAllExperimentsP60(:);
    k = char(Mutants{j});
    k1 = strrep(k,'_','-');
    title(strcat('Bulbous-',k1),'FontSize',16)
    if j == 1
        ylabel('Probability','FontWeight','bold','FontSize',14)
    end
    %ylabel('P60','FontSize',16)
    [counts,edges] = histcounts(Data,edges);
    histogram('BinEdges',edges,'BinCounts',counts./sum(counts))
    MeanL = mean(Data);
    StdL = std(Data); 
    line([MeanL  MeanL],[0 0.2],'Color','r','LineStyle','--','LineWidth',3)
    line([MeanL+StdL MeanL+StdL],[0 0.2],'Color','r','LineStyle',':','LineWidth',2)
    line([MeanL-StdL  MeanL-StdL],[0 0.2],'Color','r','LineStyle',':','LineWidth',2)
    ylim([0 0.7])
    xlabel('Number','FontSize',16)
                     
   
    %print(201+j,'-dpng',strcat('FiloNumberDistributions',char(Mutants{j}),'.png'))

    
    figure(101)
    hold on
    times = 1:60;
    %Nr Filo P60
    Means1 = mean(NbrsBulbousAllExperimentsP60);
    %Means1(TF) = nan;
    plot(times,Means1,'k-','LineWidth',4)
    ylim([0 5.5])
    if j == 1
        ylabel('Number','FontWeight','bold','FontSize',14)
    end
    k = char(Mutants{j});
    k1 = strrep(k,'_','-');
    title(strcat('Bulbous-',k1),'FontSize',16)
    xlabel('time (min)','FontWeight','bold','FontSize',12)
   % print(101+j,'-dpng',strcat('FiloNumbers',char(Mutants{j}),'.png'))

end % end over Mutants

save('./AllData.mat','Output');
