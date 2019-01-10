%% Author : Romain Fournier
% Last Modification : 10.01.2019
% Goal : Get the legendre coefficients of a Monte Carlo simulation 

%% Import the data 
% Description : Currently, the monte carlo script writes one file per
% timestep, due to the lazyness of the developper to write in parallel.
% Therefore, we must import all file to get the time-correlation function

% Tune the variables to your particular situation

% 1 General 
f=figure;
check_plot = true ; % True if we want to see the data
simulation_name="harmonic_oscillator"; % for figure names
save_folder = "harm_osc_results/"; % Folder in which we save the figures and coefficients
%2 MonteCarlo
test = false; % if true, we generate data from exponential
folder = "../MaxEnt/harm_osc/"; % Folder in which the files are 
filename_prefix = "simulation"; % Prefix (common to all files)
nb_timesteps = 31 ; % The total number of timesteps
%3 Legendre parameters
NB_GLs = [4,8,16,32,64]; % array of number of coefficients to test (to know how many of them we need)
% No need to modify the next lines 
mean_data=[];
var_data=[];
tau=[];

if(~test)
    for i=0:nb_timesteps
        filename=folder+filename_prefix+int2str(i)+".csv";
        tmp=importdata(filename);
        mean_data=[mean_data, tmp(:,2)];
        if(length(tmp(1,:))<3)
               var_data=[var_data,0];
        else
            var_data=[var_data,tmp(:,3)];
        end
        tau=[tau,(tmp(:,1))];
    end
else
    tau=linspace(0,0.5,nb_timesteps+1);
    mean_data=0.05*exp(-10*tau);
    var_data=repmat(0.001,length(tau),1);
end

if(check_plot)
    subplot(1,3,1)
    hold on
    errorbar(tau,mean_data,(var_data))
    xlabel("\tau")
    ylabel("C^i")
    title("QMC simulation")
    
    subplot(1,3,3)
    hold on
    errorbar(tau,mean_data,(var_data))
    title("Reconstruction")
end

%% Legendre transform 
% We have the function C^i(\tau), but the neural networks take its legendre
% coefficients as input. We need to get them. 
for NB_GL = NB_GLs
    % We need to integrate between -1 and 1, so we need to rescale tau
    tau_rescaled = 2*(tau-tau(1))/(tau(end)-tau(1))-1;
    nl_handler=@(l) trapz(mean_data.*legendreP(l,tau_rescaled))*(tau_rescaled(2)-tau_rescaled(1));
    nl=zeros(NB_GL,1);
    for i  = (0:(NB_GL-1))
       nl(i+1)= nl_handler(i);
    end
    if(check_plot)
        subplot(1,3,2)
        plot(nl)
        xlabel("Coefficient")
        ylabel("G_l")
        title("Legendre Coef")
    end

    %% Reconstruct the time correlation function from legendre coefficient
    % We want to make sure that we have "good" information about the original
    % function

    c_i_reconstructed = zeros(1,nb_timesteps+1);
    for i = 0 : NB_GL-1
        c_i_reconstructed = c_i_reconstructed +(2*i+1)* nl(i+1)*legendreP(i,tau_rescaled)/2;
    end
    if(check_plot)
        subplot(1,3,3)
        hold on
        plot(tau,c_i_reconstructed)
        xlabel("\tau")
        ylabel("c_i Legendre")
        title("Imaginary Time Correlation")
    end
end
%% Save the figure
if(check_plot)
    subplot(1,3,3)
    legendCell = [{"QMC"} strcat('nb coefs =',string(num2cell(NB_GLs)))];
    legend(legendCell);
    print_figure(f,save_folder+simulation_name+"_correlation",24,16);
end

