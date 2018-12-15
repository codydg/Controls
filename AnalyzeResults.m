%%
%

if ~exist('CombinedResults','var')
    uiload;
end

summary = zeros(size(CombinedResults, 1), 11);
for i = 1:size(summary,1)
    linResult = CombinedResults{i,1};
    nonLinResult = CombinedResults{i,2};
    R = CombinedResults{i,3};
    Q = CombinedResults{i,4};
    
    linResult = linResult(end-100:end-1,[2 4 6]);
    nonLinResult = nonLinResult(end-100:end-1,[2 4 6]);
    
    summary(i,1) = mean(abs(linResult(:)));
    summary(i,2) = max(abs(linResult(:)));
    summary(i,3) = mean(abs(nonLinResult(:)));
    summary(i,4) = max(abs(nonLinResult(:)));
    summary(i,5)  = R;
    summary(i,6)  = Q(1,1);
    summary(i,7)  = Q(2,2);
    summary(i,8)  = Q(3,3);
    summary(i,9)  = Q(4,4);
    summary(i,10) = Q(5,5);
    summary(i,11) = Q(6,6);
end

[~,idxLin] = sort(summary(:,1));
[~,idxNonLin] = sort(summary(:,2));

subplot(4,4,1); plot(log10(summary(idxLin,5))); title('R');
subplot(4,4,2); plot(log10(summary(idxLin,6))); title('Q1');
subplot(4,4,5); plot(log10(summary(idxLin,7))); title('Q2');
subplot(4,4,6); plot(log10(summary(idxLin,8))); title('Q3');
subplot(4,4,9); plot(log10(summary(idxLin,9))); title('Q4');
subplot(4,4,10); plot(log10(summary(idxLin,10))); title('Q5');
subplot(4,4,13); plot(log10(summary(idxLin,11))); title('Q6');

subplot(4,4,3); plot(log10(summary(idxNonLin,5))); title('R');
subplot(4,4,4); plot(log10(summary(idxNonLin,6))); title('Q1');
subplot(4,4,7); plot(log10(summary(idxNonLin,7))); title('Q2');
subplot(4,4,8); plot(log10(summary(idxNonLin,8))); title('Q3');
subplot(4,4,11); plot(log10(summary(idxNonLin,9))); title('Q4');
subplot(4,4,12); plot(log10(summary(idxNonLin,10))); title('Q5');
subplot(4,4,15); plot(log10(summary(idxNonLin,11))); title('Q6');