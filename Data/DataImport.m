close all
clear all


path = './Files/';
Mutants = {'Control';'Atg6_mutant'; 'Atg7_mutant'; 'Atg6_rescue'};

LClassificationThreshold = 8% in min shortlived have lifetime smaller than this threshold  change in lines 113 and 464


MeanNbrsFilos = zeros(length(Mutants),1);
SDNbrsFilos = zeros(length(Mutants),1);
MeanNbrsStabFilos = zeros(length(Mutants),1);
SDNbrsStabFilos = zeros(length(Mutants),1);

for j = 1:length(Mutants)
    MutantName = strrep(char(Mutants{j}),'_','-');
    %Path name and entering the path
    FullPath = strcat(path,char(Mutants{j}));
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
        %Asssign files to P40 or P60

        % save statistics for each mutant
        NrP60 = NrFiles;% Nr of P60 files
        FiloInfoP60 = nan(300,4,NrP60);%For each Filo (dim1): Start, End, Lifetime (dim2), dim3 = growth cone
        FiloInfoP60_F = nan(300,4,NrP60);
        FiloInfoP60_sF = nan(300,4,NrP60);

        LifeTimeValues = nan(500,NrP60,1);
        LifeTimeValuesStab = nan(500,NrP60,1);
        
        NbrsStabFiloAllExperimentsP60 = zeros(NrP60,60);
        NbrsFilosAllExperimentsP60 = zeros(NrP60,60);
        counterP60 = 0;
        for z = 1:length(Filenames)
            filename = char(Filenames{z});
            GCPos = strfind(filename,'GC')+2;
            GC = str2num(filename(GCPos));
            %read table
            T = readtable(filename);
            NrEntries = length(T.StartTimeStep);

            StartTimes = T.StartTimeStep;
            EndTimes = T.EndTimeStep;
            try
                LifeTimes = T.LifeTime_min_;
            catch
                LifeTimes = T.lifetime;
            end
            try
                IsBulbous = logical(T.Bulbous);
            catch
                IsBulbous = false(NrEntries,1);
            end
            %save the data

            FiloInfoP60(1:NrEntries,1,z) = StartTimes;
            FiloInfoP60(1:NrEntries,2,z) = EndTimes;
            FiloInfoP60(1:NrEntries,3,z) = LifeTimes;
            FiloInfoP60(1:NrEntries,4,z) = GC;

            % 2c. remove outliers; outliers are most likely stabilized
            % Filopodia or bulbous tips
            TF_fast = (LifeTimes < LClassificationThreshold) & ~IsBulbous;
            TF_stable = (LifeTimes >= LClassificationThreshold) & ~IsBulbous;
            %MinOutlierTime = min(LifeTimes(TF));
            LifeTimesNew = LifeTimes(TF_fast);
            LifeTimesStab = LifeTimes(TF_stable);
            
            LifeTimeValues(1:length(LifeTimesNew),z) = LifeTimesNew;
            LifeTimeValuesStab(1:length(LifeTimesStab),z) = LifeTimesStab;
            %LifeTimeValuesBulb(1:length(LifeTimesBulb),z) = LifeTimesBulb;
            
             %% Compute a number of statistics
            %Plot Filo & Stabous Numbers
            FiloInfoTmp1 = FiloInfoP60(TF_fast,:,z);
            FiloInfoTmp2 = FiloInfoP60(TF_stable,:,z);
            
            FiloInfoP60_F(1:sum(TF_fast),:,z) = FiloInfoTmp1;
            FiloInfoP60_sF(1:sum(TF_stable),:,z) = FiloInfoTmp2;
            
            NrFilos = sum(TF_fast);
            NrStab = sum(TF_stable);
            AllNrFilos = zeros(1,60);
            for counter = 1:NrFilos
                Start = FiloInfoTmp1(counter,1)+1;
                End = FiloInfoTmp1(counter,2)+1;
                AllNrFilos(Start:End) = AllNrFilos(Start:End)+1;
            end
            AllNrStab = zeros(1,60);   
            for counter = 1:NrStab
                Start = FiloInfoTmp2(counter,1)+1;
                End = FiloInfoTmp2(counter,2)+1;
                AllNrStab(Start:End)= AllNrStab(Start:End)+ 1;
            end
            figure(101+j)
            subplot(1,3,1)
            hold on
            plot(AllNrFilos,':','LineWidth',2)
            subplot(1,3,2)
            hold on
            plot(AllNrStab,':','LineWidth',2)
            %------------------
            counterP60 = counterP60 +1;
            NbrsStabFiloAllExperimentsP60(counterP60,:) = AllNrStab;
            NbrsFilosAllExperimentsP60(counterP60,:) = AllNrFilos;

        end %end over files
        % --- save data in data structure
        if z ~= 3
            disp('correct the code here')
            return
        end
        FiloInfoP60_F_tmp = [FiloInfoP60_F(:,:,1);FiloInfoP60_F(:,:,2);FiloInfoP60_F(:,:,3)];
        LGX = ~isnan(FiloInfoP60_F_tmp(:,1));
        FiloInfoP60_F_tmp2 = FiloInfoP60_F_tmp(LGX,:);
        Output.(char(Mutants(j))).P60.F.LTimes = FiloInfoP60_F_tmp2(:,3);    
        Output.(char(Mutants(j))).P60.F.StartTimes = FiloInfoP60_F_tmp2(:,1); 
        Output.(char(Mutants(j))).P60.F.EndTimes = FiloInfoP60_F_tmp2(:,2); 
        Output.(char(Mutants(j))).P60.F.GC = FiloInfoP60_F_tmp2(:,4); 
        
        FiloInfoP60_sF_tmp = [FiloInfoP60_sF(:,:,1);FiloInfoP60_sF(:,:,2);FiloInfoP60_sF(:,:,3)];
        LGX = ~isnan(FiloInfoP60_sF_tmp(:,1));
        FiloInfoP60_sF_tmp2 = FiloInfoP60_sF_tmp(LGX,:);
        Output.(char(Mutants(j))).P60.sF.LTimes = FiloInfoP60_sF_tmp2(:,3);    
        Output.(char(Mutants(j))).P60.sF.StartTimes = FiloInfoP60_sF_tmp2(:,1); 
        Output.(char(Mutants(j))).P60.sF.EndTimes = FiloInfoP60_sF_tmp2(:,2); 
        Output.(char(Mutants(j))).P60.sF.GC = FiloInfoP60_sF_tmp2(:,4); 
        %-----
    cd ..
    cd ..
    
    figure(201+j)
    hold on
    edges = (0:20)-0.5;
    for i = 1:2
        subplot(1,2,i)
         hold on
        switch i
            case 1
                Data = NbrsFilosAllExperimentsP60(:);
                title(strcat('Nr. (detect) filo.(',MutantName,')'),'FontSize',16)
                ylabel('P60','FontSize',16)
            case 2
                 title(strcat('Nr. stab. filo.(',MutantName,')'),'FontSize',16)
                 Data = NbrsStabFiloAllExperimentsP60(:);
               xlabel('number','FontSize',14)
        end
         [counts,edges] = histcounts(Data,edges);
        histogram('BinEdges',edges,'BinCounts',counts./sum(counts))
        MeanL = mean(Data);
        StdL = std(Data); 
        line([MeanL  MeanL],[0 0.2],'Color','r','LineStyle','--','LineWidth',3)
        line([MeanL+StdL MeanL+StdL],[0 0.2],'Color','r','LineStyle',':','LineWidth',2)
        line([MeanL-StdL  MeanL-StdL],[0 0.2],'Color','r','LineStyle',':','LineWidth',2)
        ylim([0 0.3])
                     
        if i == 1 
            MeanNbrsFilos(j) = MeanL;
            SDNbrsFilos(j) =StdL;
        else
           MeanNbrsStabFilos(j) = MeanL;
            SDNbrsStabFilos(j) = StdL;
        end
    end
   
    %print(201+j,'-dpng',strcat('FiloNumberDistributions',char(Mutants{j}),'.png'))

    
    figure(101+j)
    subplot(1,3,1)
    times = 1:60;
    %Nr Filo P60
    Means1 = mean(NbrsFilosAllExperimentsP60);
    TF = isoutlier(Means1);
    Means1(TF) = nan;
    plot(times,Means1,'k-','LineWidth',4)
    ylim([0 15])
    %Nr Stab. Filo P60
    Means3 = mean(NbrsStabFiloAllExperimentsP60);
    TF = isoutlier(Means3);
    Means3(TF) = nan;
    subplot(1,3,2)
    plot(times,Means3,'k-','LineWidth',4)
    ylim([0 20])

    %Ratio in P60
    Ratio1 = Means3./Means1;% correction for non-detection
    TF = isoutlier(Ratio1,'movmedian',5);
    Ratio1(TF) = nan;
    Ratio1(isnan(Ratio1)) = inf;
    TF2 = isinf(Ratio1);
    Ratio1 = Ratio1(~TF2);
    subplot(1,3,3)
    plot(times(~TF2),Ratio1,'k-','LineWidth',4)
    hold on
    subplot(1,3,3)
    Ratio1(isnan(Ratio1)) = inf;
    Ratio1 = Ratio1(~isinf(Ratio1));
    line([1 60],[mean(Ratio1) mean(Ratio1)],'LineStyle',':','Color','r','LineWidth',3)
    ylim([0 10])
    title('Ratio: stab.Filo-to-Filo')
    
    figure(101+j)
    subplot(1,3,1)
    ylabel('P60','FontWeight','bold','FontSize',14)
    title(strcat('Nr. (detectable) filo. (',MutantName,')'))
    xlabel('time (min)','FontWeight','bold','FontSize',12)
    subplot(1,3,2)
    title(strcat('Nr. stab. filo.(',MutantName,')'))
    xlabel('time (min)','FontWeight','bold','FontSize',12)
    subplot(1,3,3)
    xlabel('time (min)','FontWeight','bold','FontSize',12)
   % print(101+j,'-dpng',strcat('FiloNumbers',char(Mutants{j}),'.png'))

end % end over Mutants

save('./AllData.mat','Output');
