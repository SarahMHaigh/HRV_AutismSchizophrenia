% Analyzing ECG peaks with stats
% From preprocessed data with artefact rejection (HRV_ECGpre_toAnalysisGit.m)

clear all
close all

path = '~/Documents/ECG/corReFilt/';
Ss = {'S2' 'S3' 'S5' 'S6' 'S7' 'S8' 'S9' 'S10' 'S11', 'S12' 'S13' 'S14'... % #12
     'A1' 'A2' 'A3' 'A4' 'A7' 'A8' 'A9' 'A10' 'A14' 'A16' 'A17' 'A18' 'A19' 'A20' 'A22' 'A23' 'A24' 'A25'... % #30 'A5' 'A6' 'A11' 'A12' 'A13' 'A15' 'A21'
     'C1' 'C2' 'C3' 'C4' 'C5' 'C6' 'C9' 'C12' 'C13' 'C14' 'C15' 'C18' 'C19' 'C20' 'C21' 'C22' 'C23' 'C24' 'C25' 'C26' 'C27' 'C28'}; % 'C7' 'C10' 'C11' 'C16' 'C17'
 
Group = [zeros(12,1); ones(18,1); repmat(2,22,1)];
szEnd = 12;
autEnd = 30;

 for i = 1:length(Ss)
    newO{i} = csvread([path Ss{i} '/' Ss{i} '_corECGonset.csv']);
    new{i} = csvread([path Ss{i} '/' Ss{i} '_corECGamp.csv']);
    check = [newO{i} new{i}(1,:)'];
    
    for j = 1:length(check)-1
        if isnan(check(j,2)) || isnan(check(j+1,2))
            RR(j) = NaN;
        else RR(j) = check(j+1,1)-check(j,1);
            RRsq(j) = RR(j)^2;
        end
    end
    
    RRpeak(i) = nanmean(RR);
    RRsd(i) = nanstd(RR);
    RRrootMsq(i) = sqrt((sum(RRsq))/j);
    
 end
 
%%%%% Data checks %%%%%%%%%%
%% Load and count
for i = 1:length(Ss)
    new{i} = csvread([path Ss{i} '/' Ss{i} '_corECGamp.csv']);
    newT{i} = csvread([path Ss{i} '/' Ss{i} '_corECGtime.csv']);
    oCount(i) = length(data);
    for j = 1:length(data)
        P = data(1,j);
        Q = data(2,j);
        R = data(3,j);
        S = data(4,j);
        T = data(5,j);
        RS = R-S;
        if P>Q && R>Q && RS>50 && T>S
            new{i}(:,j) = data(:,j);
            newT{i}(:,j) = timeD(:,j);
        else new{i}(1:5,j) = NaN;
            newT{i}(1:5,j) = NaN;
        end
    end
    nCount{1,i} = Ss(i);
    nCount{2,i} = nnz(~isnan(new{i}(1,:)));
    nCount{3,i} = oCount(1,i) - nCount{2,i};
    if nCount{3,i} / nCount{2,i} > 0.5
        nCount{1,i}
    end
    clearvars j P Q R S T RS
end

%% Collating peak amplitudes and SD
for i = 1:length(Ss)
    Pmean(i,:) = nanmean(new{i}(1,:));
    Pstd(i,:) = nanstd(new{i}(1,:));
    Qmean(i,:) = nanmean(new{i}(2,:));
    Qstd(i,:) = nanstd(new{i}(2,:));
    Rmean(i,:) = nanmean(new{i}(3,:));
    Rstd(i,:) = nanstd(new{i}(3,:));
    Smean(i,:) = nanmean(new{i}(4,:));
    Sstd(i,:) = nanstd(new{i}(4,:));
    Tmean(i,:) = nanmean(new{i}(5,:));
    Tstd(i,:) = nanstd(new{i}(5,:));
end

%% Analysis of peak timings and SD of timings
for i = 1:length(Ss)
    for j = 1:length(newT{i})
        PR{i}(1,j) = newT{i}(3,j)-newT{i}(1,j);
        QR{i}(1,j) = newT{i}(3,j)-newT{i}(2,j);
        RS{i}(1,j) = newT{i}(4,j)-newT{i}(3,j);
        QS{i}(1,j) = newT{i}(4,j)-newT{i}(2,j);
        ST{i}(1,j) = newT{i}(5,j)-newT{i}(4,j);
        
    end
    
    PR_mean(i) = nanmean(PR{i}(1,:));
    PR_std(i) = nanstd(PR{i}(1,:));
    QR_mean(i) = nanmean(QR{i}(1,:));
    QR_std(i) = nanstd(QR{i}(1,:));
    RS_mean(i) = nanmean(RS{i}(1,:));
    RS_std(i) = nanstd(RS{i}(1,:));
    QS_mean(i) = nanmean(QS{i}(1,:));
    QS_std(i) = nanstd(QS{i}(1,:));
    ST_mean(i) = nanmean(ST{i}(1,:));
    ST_std(i) = nanstd(ST{i}(1,:));
    
    PmeanT(i,:) = nanmean(newT{i}(1,:));
    PstdT(i,:) = nanstd(newT{i}(1,:));
    QmeanT(i,:) = nanmean(newT{i}(2,:));
    QstdT(i,:) = nanstd(newT{i}(2,:));
    RmeanT(i,:) = nanmean(newT{i}(3,:));
    RstdT(i,:) = nanstd(newT{i}(3,:));
    SmeanT(i,:) = nanmean(newT{i}(4,:));
    SstdT(i,:) = nanstd(newT{i}(4,:));
    TmeanT(i,:) = nanmean(newT{i}(5,:));
    TstdT(i,:) = nanstd(newT{i}(5,:));
end
