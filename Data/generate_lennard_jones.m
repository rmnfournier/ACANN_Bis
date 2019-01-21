%% By Romain Fournier
% Last Edit : 21.01.19
% Objective : Generate Power Spectra and imaginary-time correlation data
% for Lennard-Jones-like fluids. 


type="training"; % training, validation or test
nb_data = 90; % Number of data to append at the end of the file
w = linspace (-50,50,2048 );

m=1;
beta=0.25;
tau= linspace(0,beta,1024)';

for i = 1:nb_data
    %% Dynamic friction Kernel. 
    % We start by defining a random dynamic friction kernel. Look at the
    % livescript to see the formula and its parameters.
    
    % Random parameters
    xi_0 = randi(601)-1; 
    a1 = randi(5e5);
    a2 = randi(1e3);
    alpha1 = randi(900)+100;
    alpha2 = randi(50)+50;
    
    omega_0=rand(1)*15;
    % Kernel definition 
    xi = @(t) xi_0*( exp(-alpha1*(f*t).^2).*(1+a1*(f*t).^4) + a2*(f*t).^4.*exp(-alpha2*(f*t).^2));

    %% Bath Power Spectrum 
    % This kernel allows to compute the power spectrum of the bath (see
    % live script)
    integrand=@(t,omega) cos(omega.*t).*xi(t);
    J = @(omega) omega.*integral(@(x)integrand(x,omega),0,10,'ArrayValued',true);

    %% "Self Energy" 
    % The solution of the problem is formulated with a "self-energy" term
    integrand_bis=@(t,omega) sin(omega.*t).*xi(t);
    im_se = -sign(w).*J(w)/m;
    re_se= sign(w)/m.*w.*integral(@(x)integrand_bis(x,w),0,10,'ArrayValued',true);

    %% Power Spectrum
    I = -1/m * coth(beta*w/2).*im_se./((w.^2-omega_0^2-re_se).^2+(im_se).^2);
    
    %% Imaginary-time Correlation function
    G = trapz(w(w>0), exp(-tau.*w(w>0)).*I(w>0),2);
    
    %% Save the data
    dlmwrite(strcat('./lennard_jones/A_',type,'.csv'),I(w>0)./sum(I(w>0)),'-append','delimiter',',')
    dlmwrite(strcat('./lennard_jones/G_',type,'.csv'),G,'-append','delimiter',',')

end