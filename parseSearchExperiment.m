function [eegEpoched, events] = parseSearchExperiment()
BAD = [];
addpath('Utilities');
Fs=500;
eegfile = 'V:\EEG\Visual Search Task\Fastovets\Session1\NeoRec_2018-08-21_13-33-50.edf';
logfile = 'V:\EEG\Visual Search Task\Fastovets\Session1\2018.08.21-13.33.54.122.log';
antfile = 'V:\EEG\Visual Search Task\Fastovets\Session1\NeoRec_2018-08-21_13-33-50_evt.edf';

events = getEventsFromLog(logfile, antfile);

if ~isempty(BAD)
    events=events(setdiff(1:length(events),BAD));
end

onsets= [events.time];
labels = [events.label];
[eeg] = ReadEDF(eegfile);
eeg = cell2mat(eeg);
eeg=eeg';
eeg = eeg(1:128,:);
ul = unique(labels);
ul = setdiff(ul,0);
winSize = 2*Fs;
step = 0.5*Fs;

sig = zeros(1,winSize);
[spcheck, f_] = get_spectrum (sig,Fs);
f2 = (f_>1)&(f_<35);
splen = length(spcheck);

for i=1:length(ul)
    lab = ul(i);
    nk=1;
    for k=1:length(onsets)
        if labels(k)==lab
            spectr_cut = zeros(128,splen);
            for ch=1:128
                onset_eeg = eeg(ch,onsets(k)+round(0.3*Fs): onsets(k+1)-1);
                onset_eeg = onset_eeg - mean(onset_eeg);
                
                cutted = cutSignal(onset_eeg, winSize, step);
                for kc = 1:size(cutted,1)
                    spectr_cut(ch,:) = spectr_cut(ch,:) + get_spectrum (cutted(kc,:),Fs);
                end
                spectr_cut(ch,:)=spectr_cut(ch,:)/kc;
            end
            spectr_cut = spectr_cut(:,f2);
            eegEpochedCut{i}(:,:,nk) = mean(spectr_cut,1);
            nk=nk+1;
            clear spectr
        end
    end
end

for i=1:7
    eegEpochedCut{i} = squeeze(eegEpochedCut{i});
end

cols=[{'red'}, {'blue'}, {'green'}, {'cyan'}, {'magenta'}, {'black'}, {'yellow'}];
 
figure;
for i=1:7
    hold on; plot(f_(f2), smooth(smooth(mean(eegEpochedCut{i},2))),cols{i});
end
end

function cutted = cutSignal(sig, winSize, step)
k=1;
for i=1:step:length(sig)-winSize+1
    cutted(k,:) = sig(i:i+winSize-1);
    k=k+1;
end
end
    