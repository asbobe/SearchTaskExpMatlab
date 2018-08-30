function [events, matevents] = getEventsFromLog(logfile, antfile)
addpath('Utilities');
[~, antheader] = ReadEDF(antfile);
Fs = 500;
antEvents = antheader.annotation.starttime(2:end);
firstShow = round(antEvents(1)*Fs);
cats{1} = 'abstract';
cats{2} = 'animals';
cats{3} = 'buildings';
cats{4} = 'face';
cats{5} = 'furniture';
cats{6} = 'kitchen';
cats{7} = 'transport';

raw = importlog(logfile);
numEvents = size(raw,1);
n=1;
messages = raw(:,4);
eventTypes = cell2mat(raw(:,2));
times = round(cell2mat(raw(:,3))/1000*Fs);
for i=1:numEvents;
    if ~isempty(strfind(eventTypes(i,:), 'SMALL_IMAGE'))
        for j=1:length(cats)
            if ~isempty(strfind(messages{i}, cats{j}))
                events(n).type = 1;
                events(n).label = j;
                events(n).time = firstShow + times(i);
                n=n+1;
            end
        end
    end
    if ~isempty(strfind(eventTypes(i,:), 'BIG_IMAGE'))
        events(n).type = 2;
        events(n).label = events(n-1).label;
        t=1;
        while ~(events(n-t+1).type==1)
            events(n).label = events(n-t).label;
            t=t+1;
        end;
        events(n).time = firstShow + times(i);
        n=n+1;
    end
    if ~isempty(strfind(eventTypes(i,:), 'PAUSE'))
        events(n).type = 0;
        events(n).label = 0;
        events(n).time = firstShow + times(i);
        n=n+1;
    end;
end

for i=1:size(events,1); matevents(i,1) = events(i).type; matevents(i,2) = events(i).label; matevents(i,3) = events(i).time; end;