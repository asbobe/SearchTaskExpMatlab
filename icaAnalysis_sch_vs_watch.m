function T = icaAnalysis_sch_vs_watch()

addpath('FastICA_2.5');
BAD = [];
addpath('Utilities');
Fs=500;
eegfile = 'E:\Databases\EEG\Visual Search Task\Alekseev\Session1\NeoRec_2018-08-23_15-41-44.edf';
logfile = 'E:\Databases\EEG\Visual Search Task\Alekseev\Session1\2018.08.23-15.42.00.825.log';
antfile = 'E:\Databases\EEG\Visual Search Task\Alekseev\Session1\NeoRec_2018-08-23_15-41-44_evt.edf';
[eeg] = ReadEDF(eegfile);

cols=[{'red'}, {'blue'}, {'green'}, {'cyan'}, {'magenta'}, {'black'}, {'yellow'}];

events = getEventsFromLog(logfile, antfile);

if ~isempty(BAD)
    events=events(setdiff(1:length(events),BAD));
end
onsets= [events.time];
labels = [events.label];

eeg = cell2mat(eeg);
eeg=eeg';

eeg = eeg(1:128,:);
k=1;
eegSearch = zeros(128,0);
for i=1:length(onsets)-1;
    if events(i).type==1;
        searchOnsets(k)=onsets(i);
        searchLabels(k) = 1;
        searchTime(k) = onsets(i+1)-onsets(i);
        eegSearch=horzcat(eegSearch,eeg(:,onsets(i):onsets(i+1)));
        k=k+1;
    end;
    if events(i).type==2;
        searchOnsets(k)=onsets(i);
        searchLabels(k) = 2;
        searchTime(k) = onsets(i+1)-onsets(i);
        eegSearch=horzcat(eegSearch,eeg(:,onsets(i):onsets(i+1)));
        k=k+1;
    end;
end;

%[dataIca, A, T] = fastica(eegSearch);
load('T_Alekseev.mat');
dataIca = T*eeg;

numComps = size(dataIca,1);
winSize = 3*Fs;
step = 0.5*Fs;

nk=ones(2);

for i=1:length(searchOnsets)
    lab = searchLabels(i);
    if searchTime(i)>winSize
        onset_eeg = dataIca(:,searchOnsets(i):searchOnsets(i)+searchTime(i)-1);
        cutted = cutSignal(onset_eeg, winSize, step);
        for ct = 1:size(cutted,3)
            cur_onset_eeg = cutted(:,:,ct);
            for ch=1:numComps                 
                sample = cur_onset_eeg(ch,:) - mean(cur_onset_eeg(ch,:));
                [spectr, f] = get_spectrum (sample,Fs);
                spectr_cut = spectr(:,(f>1)&(f<35));
                f1 = f((f>1)&(f<35));       
                eegEpoched{lab}(ch,:,nk(lab)) = mean(spectr_cut,1);
            end
            nk(lab)=nk(lab)+1;
        end
    end
end

for i=1:2
    meanEpochSpectr{i} = mean(eegEpoched{i},3);
end

for i=1:numComps
    figure; title(mat2str(i-1));
    for k=1:2  
        hold on; plot(f1, smooth(smooth(meanEpochSpectr{k}(i,:))),cols{k});
    end
end
end

function cutted = cutSignal(sig, winSize, step)
k=1;
for i=1:step:size(sig,2)-winSize+1
    cutted(:,:,k) = sig(:,i:i+winSize-1);
    k=k+1;
end
end

%Alekseev
%8   9.5-12.5
%26  7.5-10
%47  5-30
%65  10-30
%79  8.5-13 
%84  3.5-5.5
%84  6-7.5
%86  6-9
%86  9-11
%88  2.5-10
%91  8-11.5
%107 8-11.5

