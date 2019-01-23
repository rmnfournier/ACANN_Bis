%% By Romain Fournier
% Last edit : 23.01.2019
% Figure describing the dataset

%% Parameters
nb_coef=10;
beta=0.25;
tau= linspace(0,beta,1024)';
g_to_plot=[17,22,71,100];
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
subplot(2,2,3)
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
    A=importdata("figure1_A.csv");
else
    G = importdata("../Data/lennard_jones/G_validation.csv");
    A=importdata("../Data/lennard_jones/A_validation.csv");
    G=G(g_to_plot,:);
    A=A(g_to_plot,:);
end
subplot(2,2,1);
plot(tau,-G,'LineWidth',1.2);
xlabel("\tau")
ylabel("G(\tau)")

subplot(2,2,2);
plot(linspace(0,25,1024/2),A(:,1:512),'LineWidth',1.2);
xlabel("\omega")
ylabel("I(\omega)")

%% Scatter plot
if replot
    G_r=importdata("figure1_G_r.csv");
else
    G_r = importdata("../Data/lennard_jones/G_validation_reduced.csv");
end

subplot(2,2,4)
scatter3(G_r(:,1),G_r(:,2),G_r(:,3),'filled')
xlabel("\lambda_1")
ylabel("\lambda_2")
zlabel("\lambda_3")
box on
grid on

%% Save figure and data
print_figure(f,'database',8.8,6.6);
if ~replot
dlmwrite('figure1_ind_contr.csv',individual_contributions,'delimiter',',');
dlmwrite('figure1_evr.csv',evr,'delimiter',',','precision','%.15f');
dlmwrite('figure1_G.csv',G,'delimiter',',');
dlmwrite('figure1_A.csv',A,'delimiter',',');

dlmwrite('figure1_tau.csv',tau,'delimiter',',');
dlmwrite('figure1_G_r.csv',G_r,'delimiter',',');
end

