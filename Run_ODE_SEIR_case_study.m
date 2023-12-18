
%% Script for running ODE model

%Set the number of symptom severity levels 
para.n_severity = 2;

%Set the number of age classes
para.n_age_class = 1;

%Define time to run model for (in days)
maxtime = 8*365;

%Define parameter options
alpha_opts = 0:0.02:1; % Strength of symptom propagation
nu_opts = 0:0.02:1; % Baseline probaility of severe disease

%Define parameters
para.gamma=[1/5 1/7]; % Mild and severe recovery rates

para.epsilon = 1/5; % Rate of becoming infectious

beta = [0.2, 0.4]; % Mild and severe transmission rates

Q = 0.5; % Quarantining rate

%Initialise initial conditions arrays
E0 = zeros(para.n_age_class, para.n_severity);
R0 = zeros(para.n_age_class, para.n_severity);

%initialise results arrays
final_size_q_all = zeros(length(nu_opts),length(alpha_opts),2);
final_size_q_sev = zeros(length(nu_opts),length(alpha_opts),2);


for nu_itr = 1:length(nu_opts)

    nu = nu_opts(nu_itr);
    
    para.nu = [1-nu, nu];
    
    I0 = 0.01*para.nu;
    
    %Set the remaining population to susceptible
    S0 = 1-sum(E0,2)-sum(I0,2)-sum(R0,2);
    
    %Define initial conditions as a structure
    ICs = struct('S',S0,'E',E0,'I',I0,'R',R0);

    for alpha_itr = 1:length(alpha_opts)

        alpha = alpha_opts(alpha_itr);

        %Value of alpha for this run
        para.alpha = alpha*ones(para.n_severity,1);

        %% Quarantining all
        para.beta = (1-Q)*beta;

        %Run model using the parameters defined above and save into a struct array
        [Classes] = ODE_SEIR_case_study(para,ICs,maxtime);
        
        %Populate results array
        final_size_q_all(nu_itr, alpha_itr,:) = Classes.R(end,:,:);

        %% Quarantining only severe
        para.beta = [beta(1),(1-Q)*beta(2)];

        %Run model using the parameters defined above and save into a struct array
        [Classes] = ODE_SEIR_case_study(para,ICs,maxtime);
        
        %Populate results array
        final_size_q_sev(nu_itr, alpha_itr,:) = Classes.R(end,:,:);

    end

end

%Defining colormaps
blue_to_red = customcolormap([0 0.4 1], {'#EE0000','#ffffff','#1854B4'});
%blue_to_red2 = customcolormap([0 0.5 1], {'#EE0000','#ffffff','#1854B4'});

%Create figure
figure(1)
set(gcf,'units','inch','position',[0,0,10,4])
tlo1 = tiledlayout(1,2);
tlo1.Padding = 'compact';

ax1 = nexttile(tlo1);
h = pcolor(nu_opts,alpha_opts,100*sum(final_size_q_all - final_size_q_sev,3));
set(h, 'EdgeColor', 'none');
clim([-30 20])
colormap(ax1, blue_to_red)
c = colorbar;
c.Ruler.TickLabelFormat='%g%%';
hold on
contour(nu_opts,alpha_opts,100*sum(final_size_q_all - final_size_q_sev,3),[-0.001,0.001],'k','LineWidth',1.5)
hold off
title({'Increase in proportion' 'infected'})
ylabel('Baseline probability, \nu')
xlabel('Dependence on infector, \alpha')

ax2 = nexttile(tlo1);
h = pcolor(nu_opts,alpha_opts,100*(final_size_q_all(:,:,2) - final_size_q_sev(:,:,2)));
set(h, 'EdgeColor', 'none');
clim([-30 20])
colormap(ax2,blue_to_red)
c = colorbar;
c.Ruler.TickLabelFormat='%g%%';
hold on
contour(nu_opts,alpha_opts,100*(final_size_q_all(:,:,2) - final_size_q_sev(:,:,2)),[-0.001,0.001],'k','LineWidth',1.5)
hold off
title({'Increase in proportion' 'severely infected'})
ylabel('Baseline probability, \nu')
xlabel('Dependence on infector, \alpha')

