function getConfusion(targets, outputs)
targets1 = zeros(max(targets), length(targets));
outputs1 = targets1;

for i=1:length(targets);
    targets1(targets(i),i) = 1;
    outputs1(outputs(i),i) = 1;
end;

plotconfusion(targets1, outputs1);