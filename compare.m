%%
%
%#ok<*SAGROW>

clc, clear, close all;

LIN_NON_BOTH = 1;
TEST_CASE = 3;
lim = 5;

load('Results.mat');
if TEST_CASE == 1
    results = results(1:4,:);
elseif TEST_CASE == 2
    results = results(5:8,:);
elseif TEST_CASE == 3
    results = results([1 8 9],:);
end

figure('units','normalized','outerposition',[0 0 1 1]);
legends = cell(0);
for i = 1:size(results,1)
    resultLin = results{i,1};
    resultNonLin = results{i,2};
    switch LIN_NON_BOTH
        case 3
            legends{end+1} = num2str(i);
            legends{end+1} = num2str(i);
            subplot 321;plot(resultLin(1:end,1),resultLin(1:end,2) - resultLin(1:end,9),resultNonLin(1:end,1),resultNonLin(1:end,2) - resultNonLin(1:end,9),'-.'); ylabel('X'); xlabel('Time'); title('Estimation Error X'); hold on; ylim([-lim,lim]);
            subplot 323;plot(resultLin(1:end,1),resultLin(1:end,4) - resultLin(1:end,11),resultNonLin(1:end,1),resultNonLin(1:end,4) - resultNonLin(1:end,11),'-.'); ylabel('Theta 1'); xlabel('Time'); title('Estimation Error Theta 1'); hold on; ylim([-lim,lim]);
            subplot 325;plot(resultLin(1:end,1),resultLin(1:end,6) - resultLin(1:end,13),resultNonLin(1:end,1),resultNonLin(1:end,6) - resultNonLin(1:end,13),'-.'); ylabel('Theta 2'); xlabel('Time'); title('Estimation Error Theta 2'); hold on; ylim([-lim,lim]);
            subplot 322;plot(resultLin(1:end,1),resultLin(1:end,3) - resultLin(1:end,10),resultNonLin(1:end,1),resultNonLin(1:end,3) - resultNonLin(1:end,10),'-.'); ylabel('d/dt X'); xlabel('Time'); title('Estimation Error d/dt X'); hold on; ylim([-lim,lim]);
            subplot 324;plot(resultLin(1:end,1),resultLin(1:end,5) - resultLin(1:end,12),resultNonLin(1:end,1),resultNonLin(1:end,5) - resultNonLin(1:end,12),'-.'); ylabel('d/dt Theta 1'); xlabel('Time'); title('Estimation Error d/dt Theta 1'); hold on; ylim([-lim,lim]);
            subplot 326;plot(resultLin(1:end,1),resultLin(1:end,7) - resultLin(1:end,14),resultNonLin(1:end,1),resultNonLin(1:end,7) - resultNonLin(1:end,14),'-.'); ylabel('d/dt Theta 2'); xlabel('Time'); title('Estimation Error d/dt Theta 2'); hold on; ylim([-lim,lim]);
        case 2
            legends{end+1} = num2str(i);
            subplot 321;plot(resultNonLin(1:end,1),resultNonLin(1:end,2) - resultNonLin(1:end,9)); ylabel('X'); xlabel('Time'); title('Estimation Error X'); hold on; ylim([-lim,lim]);
            subplot 323;plot(resultNonLin(1:end,1),resultNonLin(1:end,4) - resultNonLin(1:end,11)); ylabel('Theta 1'); xlabel('Time'); title('Estimation Error Theta 1'); hold on; ylim([-lim,lim]);
            subplot 325;plot(resultNonLin(1:end,1),resultNonLin(1:end,6) - resultNonLin(1:end,13)); ylabel('Theta 2'); xlabel('Time'); title('Estimation Error Theta 2'); hold on; ylim([-lim,lim]);
            subplot 322;plot(resultNonLin(1:end,1),resultNonLin(1:end,3) - resultNonLin(1:end,10)); ylabel('d/dt X'); xlabel('Time'); title('Estimation Error d/dt X'); hold on; ylim([-lim,lim]);
            subplot 324;plot(resultNonLin(1:end,1),resultNonLin(1:end,5) - resultNonLin(1:end,12)); ylabel('d/dt Theta 1'); xlabel('Time'); title('Estimation Error d/dt Theta 1'); hold on; ylim([-lim,lim]);
            subplot 326;plot(resultNonLin(1:end,1),resultNonLin(1:end,7) - resultNonLin(1:end,14)); ylabel('d/dt Theta 2'); xlabel('Time'); title('Estimation Error d/dt Theta 2'); hold on; ylim([-lim,lim]);
        case 1
            legends{end+1} = num2str(i);
            subplot 321;plot(resultLin(1:end,1),resultLin(1:end,2) - resultLin(1:end,9)); ylabel('X'); xlabel('Time'); title('Estimation Error X'); hold on; ylim([-lim,lim]);
            subplot 323;plot(resultLin(1:end,1),resultLin(1:end,4) - resultLin(1:end,11)); ylabel('Theta 1'); xlabel('Time'); title('Estimation Error Theta 1'); hold on; ylim([-lim,lim]);
            subplot 325;plot(resultLin(1:end,1),resultLin(1:end,6) - resultLin(1:end,13)); ylabel('Theta 2'); xlabel('Time'); title('Estimation Error Theta 2'); hold on; ylim([-lim,lim]);
            subplot 322;plot(resultLin(1:end,1),resultLin(1:end,3) - resultLin(1:end,10)); ylabel('d/dt X'); xlabel('Time'); title('Estimation Error d/dt X'); hold on; ylim([-lim,lim]);
            subplot 324;plot(resultLin(1:end,1),resultLin(1:end,5) - resultLin(1:end,12)); ylabel('d/dt Theta 1'); xlabel('Time'); title('Estimation Error d/dt Theta 1'); hold on; ylim([-lim,lim]);
            subplot 326;plot(resultLin(1:end,1),resultLin(1:end,7) - resultLin(1:end,14)); ylabel('d/dt Theta 2'); xlabel('Time'); title('Estimation Error d/dt Theta 2'); hold on; ylim([-lim,lim]);
    end
end
subplot 321;legend(legends);
subplot 322;legend(legends);
subplot 323;legend(legends);
subplot 324;legend(legends);
subplot 325;legend(legends);
subplot 326;legend(legends);