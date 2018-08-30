function featureWatcher()

addpath('FastICA_2.5');
BAD = [];
chosenComp = [0,3,5,8,10,12,14,17,18,23,25,31,33,34,37,38,38,39,54,57,59,60,62,63,79,81,81,83,86,91,95,99,104,105,110,112,119]+1;
chosen_diap = [[2,5];[1,15];[15,30];[1,5];[1,5];[2,12];[1,5];[1,4];[1,11];[1,9];[1,3.5];[1,9];[13,25];[1,7];[1,9];[7,10.5];[8,12];[2,12];[7.5,10.5];[1,6];[8,11.5];[15,30];[5,17];[8,11.5];[17,20];[7,12];[20,23.5];[6,12];[8,10];[7,11];[2,6];[7.5,10];[7,10.5];[3,8.5];[7.5,10];[6,9];[5,10]];
addpath('Utilities');
Fs=500;
eegfile = 'E:\Databases\EEG\Visual Search Task\Fastovets\Session1\NeoRec_2018-08-21_13-33-50.edf';
logfile = 'E:\Databases\EEG\Visual Search Task\Fastovets\Session1\2018.08.21-13.33.54.122.log';
antfile = 'E:\Databases\EEG\Visual Search Task\Fastovets\Session1\NeoRec_2018-08-21_13-33-50_evt.edf';
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

load('T_Fastovets.mat');
eeg = eeg(1:128,:);
dataIca = T*eeg;
winSize = 3*Fs;
step = 1*Fs;
for n=1:length(chosenComp)
    data = dataIca(chosenComp(n),:);
    graphPt = 1;
    figure; title([mat2str(chosenComp(n)-1) '; diap: ' mat2str(chosen_diap(n,:))]);
    for i=1:length(labels)-1
        lab = labels(i);
        if ~(lab==0) && onsets(i+1)-onsets(i)-0.3*Fs>winSize
            onset_data = data(onsets(i)+round(0.3*Fs): onsets(i+1)-1);
            onset_data = onset_data-mean(onset_data);
            cutted = cutSignal(onset_data, winSize, step);
            for nc = 1:size(cutted,1);
                [spectr, f] = get_spectrum (cutted(nc,:),Fs);
                specVal(nc) = mean(spectr(:,(f>chosen_diap(n,1))&(f<chosen_diap(n,2))));
            end
            hold on; plot(graphPt:graphPt+nc-1, specVal, cols{lab});
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