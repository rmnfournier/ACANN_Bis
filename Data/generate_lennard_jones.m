%% By Romain Fournier
% Last Edit : 21.01.19
% Objective : Generate Power Spectra and imaginary-time correlation data
% for Lennard-Jones-like fluids. 


type="validation"; % training, validation or test
nb_data = 100; % Number of data to append at the end of the file
w = linspace (0,50,1024 )';
good = 0;
bad=0;
m=1;
beta=0.25;
tau= linspace(0,beta/2,1024);
my_waitbar = waitbar(0,"generating data");

i=0;
for i =1:nb_data
    %% Dynamic friction Kernel. 
    % We start by defining a random dynamic friction kernel. Look at the
    % livescript to see the formula and its parameters.
    waitbar((i-1)/nb_data,my_waitbar,"generating data("+int2str(i)+"/"+int2str(nb_data)+")");
    % Random parameters
    xi_0 = randi(200)+50; 
    a1 = randi(1e5)+1.35e5;
    a2 = randi(1.5e2)+150;
    alpha1 = randi(500)+800;
    alpha2 = randi(30)+70;
    f=0.15+rand(1);
    omega_0=rand(1)*10+10;
    omega_tilde=sqrt(omega_0)-xi_0;

    % Kernel definition 
    xi = @(t) xi_0*( exp(-alpha1*(f*t).^2).*(1+a1*(f*t).^4) + a2*(f*t).^4.*exp(-alpha2*(f*t).^2));

    %% Bath Power Spectrum 
    % This kernel allows to compute the power spectrum of the bath (see
    % live script)
    integrand=@(t,omega) cos(omega.*t).*xi(t);
    J = @(omega) omega.*integral(@(x)integrand(x,omega),0,inf,'ArrayValued',true);
    %% Self Energy
    gamma = @(omega) integral(@(t) xi(t).*exp(1i*omega.*t),0,inf,'ArrayValued',true);
    g = gamma(w);
    I_handler = @(omega) sign(omega).*1/m^2.*coth(beta*omega/2).*omega.*real(g)./((omega.^2-omega_tilde^2-sign(omega).*omega/m.*imag(g)).^2+(omega/m.*real(g)).^2);
    I=I_handler(w);
    %% Imaginary-time Correlation function
   
    G = integral(@(omega) exp(-tau.*omega).*I_handler(omega),0,50,'ArrayValued',true);
    
    %% Save the data
    if(sum(isnan(I(w>0)./sum(I(w>0))))==0 && all(real(g)>=0) && all(G(2:end)>0))
        good = good+1;
        dlmwrite(strcat('./lennard_jones/A_',type,'.csv'),I'./sum(I),'-append','delimiter',',')
        dlmwrite(strcat('./lennard_jones/G_',type,'.csv'),G','-append','delimiter',',')
        dlmwrite(strcat('./lennard_jones/parameters_',type,'.csv'),[xi_0,a1,a2,alpha1,alpha2,f,omega_0],'-append','delimiter',',')
    else
        dlmwrite(strcat('./lennard_jones/err_A_',type,'.csv'),I'./sum(I),'-append','delimiter',',')
        dlmwrite(strcat('./lennard_jones/err_G_',type,'.csv'),G','-append','delimiter',',')
        dlmwrite(strcat('./lennard_jones/err_parameters_',type,'.csv'),[xi_0,a1,a2,alpha1,alpha2,f,omega_0],'-append','delimiter',',')
        bad=bad+1;
        %nb_data=nb_data+1;
    end
end
good/(good+bad)
close (my_waitbar);