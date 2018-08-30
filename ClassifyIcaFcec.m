function ClassifyIcaFcec()

addpath('FastICA_2.5');
BAD = [73];
addpath('Utilities');
Fs=500;
eegfile = 'V:/EEG/Video Observing Task/A.Alekseev/NeoRec_2018-07-13_12-06-56.edf';
labels = dlmread('V:\EEG\Video Observing Task\labels.txt');
onsets = dlmread('V:/EEG\Video Observing Task\onsets.txt');
[eeg] = ReadEDF(eegfile);

cols=[{'red'}, {'blue'}, {'green'}, {'cyan'}, {'magenta'}, {'black'}, {'yellow'}];

if ~isempty(BAD)
    labels=labels(setdiff(1:length(labels),BAD));
    onsets=onsets(setdiff(1:length(onsets),BAD));
end

onsets(end+1) = onsets(end)+500*10;
eeg = cell2mat(eeg);
eeg=eeg';

load('T_Alekseev.mat');
%[dataIca, A, T] = fastica(eeg);
dataIca = T*eeg;
ul = unique(labels);

comps = [8, 15, 33, 57, 57, 42, 56, 66, 32, 74, 72, 66, 65, 65, 19, 29, 35, 46, 53, 88, 54, 93]+1;
diaps = [[1 7]; [1 7]; [8 10]; [6.5 10]; [10.5 12]; [10 11]; [10.5 11]; [10.5 11]; [10.5 12]; [10.5 12]; [11 11.5];
    [6 11]; [7 9]; [10.5 11.5]; [10 11.5]; [11 12]; [9 11]; [10.5 12]; [9 10.5]; [11 12]; [5.5 7.5]; [10 11.3]];

for i=1:length(ul)
    lab = ul(i);
    curOnsets = onsets(labels==lab);
    for k=1:length(curOnsets)
        fvec=[];
        for n=1:length(comps)
            ch = comps(n);
            onset_eeg = dataIca(ch,curOnsets(k): curOnsets(k)+10*Fs-1);
            onset_eeg = onset_eeg - mean(onset_eeg);
            [cspectr, f] = get_spectrum (onset_eeg,Fs);
            cspectr = mean(cspectr(f>diaps(n,1)&(f<diaps(n,2))));
            fvec = horzcat(fvec, cspectr);
        end
        dataEpoched{i}(:,k) = fvec;
        clear fvec
    end
end

figure;
for i=1:7
    hold on; plot(mean(dataEpoched{i},2),cols{i});
end
