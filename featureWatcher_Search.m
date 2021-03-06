function featureWatcher_Search()

addpath('FastICA_2.5');
BAD = [];
chosenComp = [6 14 23 29 35 52 74 80 89 103]+1;
chosen_diap = [[1,15];[1,12];[20,28];[7,11];[4.5,13];[7,9];[8,13];[7,11];[4,9];[8,11]];
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

load('T_Alekseev_Search.mat');
eeg = eeg(1:128,:);
k=1;
for i=1:length(onsets)-1;
    if events(i).type==2;
        searchOnsets(k)=onsets(i);
        searchLabels(k) = labels(i);
        searchTime(k) = onsets(i+1)-onsets(i);
        k=k+1;
    end;
end;
dataIca = T*eeg;
winSize = 3*Fs;
step = 0.5*Fs;
for n=1:length(chosenComp)
    data = dataIca(chosenComp(n),:);
    graphPt = 1;
    figure; title([mat2str(chosenComp(n)-1) '; diap: ' mat2str(chosen_diap(n,:))]);
    for i=1:length(searchOnsets)-1
        lab = searchLabels(i);
        if searchTime(i)>winSize
            onset_data = data(searchOnsets(i):searchOnsets(i)+searchTime(i)-1);
            %onset_data = onset_data-mean(onset_data);
            cutted = cutSignal(onset_data, winSize, step);
            for nc = 1:size(cutted,1);
                [spectr, f] = get_spectrum (cutted(nc,:),Fs);
                specVal(nc) = mean(spectr(:,(f>chosen_diap(n,1))&(f<chosen_diap(n,2))));
            end
            hold on; plot(graphPt:graphPt+nc-1, specVal, 's', 'MarkerSize',5, 'MarkerFaceColor', cols{lab},'MarkerEdgeColor', cols{lab});
            graphPt=graphPt+nc;
            clear specVal;
        end
    end
end
end

function cutted = cutSignal(sig, winSize, step)
k=1;
for i=1:step:length(sig)-winSize+1
    cutted(k,:) = sig(i:i+winSize-1);
    k=k+1;
end
end
