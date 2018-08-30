function T = icaAnalysis_general()

addpath('FastICA_2.5');
BAD = [];
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

eeg = eeg(1:128,:);
%[dataIca, A, T] = fastica(eeg(:,onsets(1):end));
load('T_Fastovets.mat');
dataIca = T*eeg;

ul = unique(labels);
ul = setdiff(ul,0);
numComps = size(dataIca,1);
winSize = 3*Fs;
step = 0.5*Fs;

for i=1:length(ul)
    lab = ul(i);
    nk=1;
    for k=1:length(onsets)-1
        if labels(k)==lab && onsets(k+1)-onsets(k)-0.3*Fs>winSize
            onset_eeg = dataIca(:,onsets(k)+round(0.3*Fs): onsets(k+1)-1);
            cutted = cutSignal(onset_eeg, winSize, step);
            for ct = 1:size(cutted,3)
                cur_onset_eeg = cutted(:,:,ct);
                for ch=1:numComps                         
                    sample = cur_onset_eeg(ch,:) - mean(cur_onset_eeg(ch,:));
                    [spectr, f] = get_spectrum (sample,Fs);
                    spectr_cut = spectr(:,(f>1)&(f<35));
                    f1 = f((f>1)&(f<35));       
                    eegEpoched{i}(ch,:,nk) = mean(spectr_cut,1);
                end       
                nk=nk+1;
            end
        end
    end
end
for i=1:7
    meanEpochSpectr{i} = mean(eegEpoched{i},3);
end

for i=1:numComps
    figure; title(mat2str(i-1));
    for k=1:7  
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

%Fastovets Session1
%38  7-10.5   4
%38  8-12     3
%54  7.5-10.5 4
%57  1-6      5,6,7
%59  8-11.5   2,4
%63  8-11.5   4
%--79  17-20    5
%81  7-12     (up-down): (4,2)-(3,5)-1-(6,7)  attention vs search
%81  20-23.5  4 attention vs search
%83  6-12     4 => 1
%86  8-10     4
%--91  7-11     4 => (1,7)
%95  2-6      3 => 6
%--99  7.5-10   3,4
%--104 7-10.5   3 => (1,7)
%105 3-8.5    5
%110 7.5-10   3
%112 6-9      4
%119 5-10     4