%% Author : Romain Fournier
% Last Modification : 10.01.2019
% Goal : Get the legendre coefficients of a Monte Carlo simulation 

%% Import the data 
% Description : Currently, the monte carlo script writes one file per
% timestep, due to the lazyness of the developper to write in parallel.
% Therefore, we must import all file to get the time-correlation function

% Tune the variables to your particular situation
folder = "../MaxEnt/harm_osc/"; % Folder in which the files are 
filename_prefix = "simulation"; % Prefix (common to all files)
nb_timesteps = 31 ; % The total number of timesteps
check_plot = true ; % True if we want to see the data
simulation_name="harmonic_oscillator"; % for figure names
save_folder = "harm_osc_results/"; % Folder in which we save the figures

% No need to modify the next lines 
mean_data=[];
var_data=[];
tau=[];
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

if(check_plot)
    f=figure;
    errorbar(tau,mean_data,(var_data))
    xlabel("\tau")
    ylabel("C^i(\tau)");
    print_figure(f,save_folder+simulation_name+"_correlation",12,6);
end

%% 