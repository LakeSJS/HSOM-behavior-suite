function process_individual(obj)
%load(subjectName)
%close all
%% Generate SNR over trials plot
for trialNo = 1:length(obj.Trials)
    SNR(trialNo) = obj.Trials(trialNo).Parameters(10);
    Code(trialNo) = obj.Trials(trialNo).Code;
end
figure
plot(SNR)
xlabel('trial #')
ylabel('SNR')
% %% Obtain psychometric curve from proportion within condition
% SNRlevels = unique(SNR);
% for SNRno = 1:length(SNRlevels)
%     numHits = sum((SNR == SNRlevels(SNRno)) .* (Code == 1));
%     hitRate(SNRno) = numHits / sum(SNR == SNRlevels(SNRno));
% end
% figure
% plot(SNRlevels,hitRate)
%% Obtain psychometric curve from binary logistic regression
hitInds = find(Code == 1); missInds = find(Code == 2);
correctRejectInds = find(Code == 3); falsePositiveInds = find(Code == 4);
%fsigm = @(param,xval) param(1)+(param(2)-param(1))./(1+10.^((param(3)-xval)*param(4)));
S2inds = find((Code == 1) + (Code == 2));
isHit = Code == 1;
% fit binary logistic regression
x = SNR(S2inds);
y = isHit(S2inds);

mdl = fitglm(x, y, "Distribution", "binomial");
xplot = linspace(min(SNR),max(SNR),1000)'; % test data
yplot = predict(mdl,xplot);
[~,x50ind] = min(abs(yplot - 0.5));
x50 = xplot(x50ind);

%[f,ax] = makeStandardFigure(1,[1 1]);
figure
hold on
plot(x/10,y,'o')
plot(xplot/10,yplot);
xlabel('% SNR')
%printFigure('Psychometric_Curve_LJS_021223_Randomized')
%% generate confusion matrix
targets = Code <= 2; %% S2 was played if Code indicates hit or miss
outputs = (Code == 1) + (Code == 4); %% S2 spout was selected if Code indicated hit or false alarm
figure
plotconfusion(targets,outputs)
%% Generate ROC curve
% find relationship between stimulus intensity and false positive
isHit = Code == 1; hitInds = find(isHit); 
isMiss = Code == 2; missInds = find(isMiss);
isCorrectReject = Code == 3; correctRejectInds = find(isCorrectReject);
isFalsePositive = Code == 4; falsePositiveInds = find(isFalsePositive);
noiseInds = [correctRejectInds falsePositiveInds];
predictions = isHit + isFalsePositive;
trueClass = isHit + isMiss;
[X,Y] = perfcurve(predictions,trueClass,1);
figure
plot(X,Y) % not enough data, try different method
%% Generate ROC curve method 2
figure
plot(SNR(noiseInds),isFalsePositive(noiseInds),'o')

end