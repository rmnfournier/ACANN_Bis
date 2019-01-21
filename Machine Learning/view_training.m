% 2243 -> [1]
% 3334  -> [2]
% 2247  -> [1 1]
% 3344 -> [2 2]
close all
parameters = sort([149768,792584,332304,2112528,6321168]);%[2243,3334,5516,9880,18608,36064

datas=[100000];
f=figure(1);
 hold on
colors=lines(length(parameters));
d=0;
for data = datas
    d=d+1;
    c=1;
    subplot(1,length(datas),d)
    hold on
    p_txt=[];
    for p = parameters
        %% Visualize training
        if(exist("MAE_validation_"+int2str(data)+"data_"+int2str(p)+"parameters.csv"))
            p_txt=[p_txt p p ];
            val=importdata("MAE_validation_"+int2str(data)+"data_"+int2str(p)+"parameters.csv");
            tr=importdata("MAE_training_"+int2str(data)+"data_"+int2str(p)+"parameters.csv");

            %% training and validation loss
            loglog(tr(:,1),movmean(tr(:,2),100),'-','color',colors(c,:))
            loglog(val(:,1),movmean(val(:,2),100),'--','color',colors(c,:))

            c=c+1;
        end
    end
         p_txt=sort(p_txt);
        legendCell = strcat('N=',string(num2cell(p_txt)));

        legend(legendCell,'Location','southwest')
        title("loss for "+int2str(data)+" data");
        xlabel("step");
        ylabel(" Loss");
        set(gca,'yscale','log')
        set(gca,'xscale','log')
end
print_figure(f,"recap_1000data",30.6,16.3);
