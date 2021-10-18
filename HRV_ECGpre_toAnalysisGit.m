%% ECG analysis
% Preprocessing EEG data (ECG channel = 134)
% Includes artefact rejection
% Timing of Q-Q (trigger)
% High and Low frequency analysis across whole scan

clear all
close all

eeglab % excluded participants where over half of ECG epochs will be rejected
Ss = {'S2' 'S3' 'S5' 'S6' 'S7' 'S8' 'S9' 'S10' 'S11', 'S12' 'S13' 'S14'... % #12
     'A1' 'A2' 'A3' 'A4' 'A7' 'A8' 'A9' 'A10' 'A14' 'A16' 'A17' 'A18' 'A19' 'A20' 'A22' 'A23' 'A24' 'A25'... % #30 'A5' 'A6' 'A11' 'A12' 'A13' 'A15' 'A21'
     'C1' 'C2' 'C3' 'C4' 'C5' 'C6' 'C9' 'C12' 'C13' 'C14' 'C15' 'C18' 'C19' 'C20' 'C21' 'C22' 'C23' 'C24' 'C25' 'C26' 'C27' 'C28'}; % 'C7' 'C10' 'C11' 'C16' 'C17'

path = '~/Documents/'; % where .bdf files are stored
eegpath = '~/Documents/MATLAB/'; % where EEGLAB is housed

for j = 1:length(Ss)
    pathup = [path Ss{j} '/SimpleTone/'];
    EEG = pop_biosig([pathup Ss{j} '_SimpleTone.bdf'], 'ref',[129 130] ,'refoptions',{'keepref' 'off'});
    EEG.setname=[Ss(j) '_SimpleTone'];
    EEG = eeg_checkset( EEG );
    EEG=pop_chanedit(EEG, 'lookup',[ eegpath 'eeglab14_1_2b/plugins/dipfit2.3/NARSAD_cap_try.ced']);
    EEG = eeg_checkset( EEG );
    EEG = pop_eegfiltnew(EEG, 0.04,100,42240,0,[],0);
    EEG = eeg_checkset( EEG );
    EEG = pop_fmrib_qrsdetect(EEG,134,'14','no');
    EEG = pop_saveset( EEG, 'filename',['ECG/epochedNoFilt/' Ss{j} '_ECGtrig.set'],'filepath',path);
    EEG = pop_loadset('filename',['ECG/epochedNoFilt/' Ss{j} '_ECGtrig.set'],'filepath',path);
    
    %% to save onset trigger timings
    event = squeeze(struct2cell(EEG.event))';
    onset = {;};

    for i = 1:length(EEG.event)
        if strmatch('14',event(i,1),'exact')
            a(1,:) = event(i,1:3);
            sz = size(onset,1);
            if sz < 1
                if event{i+1,3}~=1
                    onset(sz+1,:) = event(i,:);
                end
            elseif onset{sz,3}~=a{1,3}
                onset(sz+1,:) = event(i,:);
            end
        end
    end
    onsetTime = cell2mat(onset(:,2));
    clearvars event onset a sz
    
    x = EEG.data(134,:);    % sampled data
    n = length(x);          % number of samples
    fs = 512;               % sampling frequency
    dt = 1/fs;              % time increment per sample
    t = (0:n-1)/fs;         % time range for data
    f = (0:n-1)*(fs/n);     % frequency range
    y = fft(x);             % discrete fourier transform (DFT) of data
    power = abs(y).^2/n;    % power of the DFT
    upperHF = 420;
    mid = 154;
    lowerLF = 42;
    figure;plot(f(lowerLF:mid),power(lowerLF:mid))
    figure;plot(f(round(mid):round(upperHF)),power(round(mid):round(upperHF)))
    xlabel('Frequency')
    ylabel('Power')
    Freq(j,1:2) = [mean(power(lowerLF:mid)) mean(power(mid:upperHF))];

    EEG = pop_epoch( EEG, {'14'}, [-0.2 0.5], 'epochinfo', 'yes');
    EEG = pop_rmbase( EEG, []);
    base = round(.1992188*512); % exact baseline in  samples
    
    ecgChan = EEG.data(134,:,:); % this number should correspond to the channel that has the ECG trace in it
    ecgChan = squeeze(ecgChan);
    
    %% find peaks, reject noisy epochs, save rest
    for i = 1:size(ecgChan,2)
        [Y, X] = max(ecgChan(75:125,i));
        ampR(i) = Y;
        timeR(i) = X+75;
        [Y1, X1] = min(ecgChan(timeR(i)-20:timeR(i),i));
        ampQ(i) = Y1;
        timeQ(i) = X1+timeR(i)-20;
        [Y, X] = min(ecgChan(120:140,i));
        ampS(i) = Y;
        timeS(i) = X+120;
        [Y, X] = max(ecgChan(150:350,i));
        ampT(i) = Y;
        timeT(i) = X+150;
        [Y, X] = max(ecgChan(20:75,i));
        ampP(i) = Y;
        timeP(i) = X+20;
    
        Q = ampQ(i);
        P = ampP(i);
        R = ampR(i);
        S = ampS(i);
        T = ampT(i);
        RS = R-S;
        beg = mean(ecgChan(1:20,i));
        ending = mean(ecgChan(338:358,i));
        diff = abs(beg-ending);
        if P>Q && R>Q && RS>50 && T>S && diff<200
            new(:,i) = ecgChan(:,i);
            nAmpQ(i) = ampQ(i);
            nAmpP(i) = ampP(i);
            nAmpR(i) = ampR(i);
            nAmpS(i) = ampS(i);
            nAmpT(i) = ampT(i);
            nTimeQ(i) = timeQ(i);
            nTimeP(i) = timeP(i);
            nTimeR(i) = timeR(i);
            nTimeS(i) = timeS(i);
            nTimeT(i) = timeT(i);
 
        else new(:,i) = NaN(EEG.pnts,1);
            nAmpQ(i) = NaN;
            nAmpP(i) = NaN;
            nAmpR(i) = NaN;
            nAmpS(i) = NaN;
            nAmpT(i) = NaN;
            nTimeQ(i) = NaN;
            nTimeP(i) = NaN;
            nTimeR(i) = NaN;
            nTimeS(i) = NaN;
            nTimeT(i) = NaN;
        end
    end
    
    clearvars P Q R S T RS ampQ ampP ampR ampS ampT timeQ timeP timeR timeS timeT beg ending diff

    figure;plot(nanmean(new,2));
    figure;plot(new);
    figure;plot(new(:,4));

    mkdir([path 'ECG/corReFilt/' Ss{j}]);
    amplitude = [nAmpP; nAmpQ; nAmpR; nAmpS; nAmpT];
    time = [nTimeP; nTimeQ; nTimeR; nTimeS; nTimeT];
    xlswrite([path 'ECG/corReFilt/' Ss{j} '/' Ss{j} '_corECGamp.xlsx'], amplitude);%'precision',6);
    xlswrite([path 'ECG/corReFilt/' Ss{j} '/' Ss{j} '_corECGtime.xlsx'], time);%'precision',6);
    xlswrite([path 'ECG/corReFilt/' Ss{j} '/' Ss{j} '_corECGecg.xlsx'], ecgChan);%'precision',6);
    xlswrite([path 'ECG/corReFilt/' Ss{j} '/' Ss{j} '_corECGnew.xlsx'], new);%'precision',6);
    xlswrite([path 'ECG/corReFilt/' Ss{j} '/' Ss{j} '_corECGonset.xlsx'], onsetTime);%,'precision',6);
    xlswrite([path 'ECG/corReFilt/' Ss{j} '/' Ss{j} '_corECGfreq.xlsx'], power);%'precision',6);
    
    clearvars nAmpP nAmpQ nAmpR nAmpS nAmpT nTimeP nTimeQ nTimeR nTimeS nTimeT amplitude time ecg new onsetTime PowerFreq
end
xlswrite([path 'ECG/corReFilt/All_corECGpower04.xlsx'], Freq);