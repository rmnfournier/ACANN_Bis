%% By Romain Fournier
% Last edit : 23.01.2019
% Figure describing the dataset

%% Parameters
nb_coef=10;
beta=0.25;
tau= linspace(0,beta,1024)';
g_to_plot=5;
f=figure;
replot=true;

%% Explain variance ratio
if replot
    evr=importdata("figure1_evr.csv");
else
    evr = importdata("../Data/lennard_jones/evr.csv");
end
individual_contributions = evr;
for i=2:length(evr)
    individual_contributions(i)=individual_contributions(i)-sum(individual_contributions(1:(i-1)));
end
subplot(1,2,1)
bar(individual_contributions(1:nb_coef));
hold on
%stairs(evr(1:nb_coef),'LineWidth',1.4)
xlabel("\lambda")
ylabel("EVR")
set(gca,'YScale','log');

%legend('Individual Contributions','Cumulative EVR','Location','northoutside')

%% Data example
if replot
    G=importdata("figure1_G.csv");
else
    G = importdata("../Data/lennard_jones/G_validation.csv");
    G=G(g_to_plot,:);
end
subplot(1,2,2);
plot(tau,-G,'LineWidth',1.2);
xlabel("\tau")
ylabel("G(\tau)")

%% Save figure and data
print_figure(f,'database',6.4,3.2);
dlmwrite('figure1_ind_contr.csv',individual_contributions,'delimiter',',');
dlmwrite('figure1_evr.csv',evr,'delimiter',',','precision','%.15f');
dlmwrite('figure1_G.csv',G,'delimiter',',');
dlmwrite('figure1_tau.csv',tau,'delimiter',',');


