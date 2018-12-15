%%
%

clc, clear, close all;

GRAPHICAL_PLOT = false;
END_PLOT = false;

initState = [ ...
    0
    0
    0
    0
    pi/16
    -pi/16
    ];

arr = [1 10 100];
i = 1;
CombinedResults = cell(0);
tic;
for R = arr
    for Q1 = arr
        for Q2 = arr
            for Q3 = arr
                for Q4 = arr
                    for Q5 = arr
                        for Q6 = arr
                            Q = diag([Q1,Q2,Q3,Q4,Q5,Q6]);
                            RunSimulation;
                            CombinedResults(end+1,:) = {controlledResultLin,controlledResultNonLin,R,Q}; %#ok<SAGROW>
                            clearvars -except R Q1 Q2 Q3 Q4 Q5 arr ...
                                initState GRAPHICAL_PLOT ...
                                END_PLOT CombinedResults;
                        end
                    end
                end
            end
        end
    end
end
toc;
[fn,pn] = uiputfile('*.mat', 'Store Results Where?');
save(fullfile(pn,fn),'CombinedResults','-v7.3');