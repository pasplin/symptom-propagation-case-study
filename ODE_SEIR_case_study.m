%{ 
ODE_SEIR_model(para,ICs,maxtime)

A function for running the ODE model for a given set of parameters

Input:
- para: A struct containing the parameters for the run
- ICs: A struct containing the initial conditions for the run
- maxtime: The time span of the simulation in days

Output:
- Classes: A struct containing the number of people in S,E,I,R at each
timestep
%}
function [Classes] = ODE_SEIR_case_study(para,ICs,maxtime)

%Run ODE using ODE45
opts = odeset('RelTol',1e-5);
tspan = 0:1:maxtime;
init_conds = [ICs.S ICs.E ICs.I ICs.R];
[t, pop] = ode45(@diff_SEIR_model, tspan , init_conds , opts, para);

%Rename these parameters to simplify code
n = 1;%para.n_age_class;
m = 2;%para.n_severity;

%Convert output to structure
Classes = struct('t',t);
Classes.S = pop(:,1:n);

%Initialise arrays
Classes.E = zeros(maxtime+1,n, m);
Classes.I = zeros(maxtime+1,n, m);
Classes.R = zeros(maxtime+1,n, m);

%Convert 2D pop array to output arrays
for ii = 1:m
    Classes.E(:,:, ii) = pop(:,ii*n+1:(ii+1)*n);
    Classes.I(:,:, ii) = pop(:,(m+ii)*n+1:(m+ii+1)*n);
    Classes.R(:,:, ii) = pop(:,(2*m+ii)*n+1:(2*m+ii+1)*n);
end

%Differential equation
function dPop = diff_SEIR_model(~,pop,para)

%Rename these parameters to simplify code
n = 1;%para.n_age_class;
m = 2;%para.n_severity;

%Convert 1D pop array to separate arrays for each class 
S=pop(1:n);

E = zeros(n, m);
I = zeros(n, m);
R = zeros(n, m);
for ii = 1:m
    E(:, ii) = pop(ii*n+1:(ii+1)*n);
    I(:, ii) = pop((m+ii)*n+1:(m+ii+1)*n);
    R(:, ii) = pop((2*m+ii)*n+1:(2*m+ii+1)*n);
end

%Initialise differential equation arrays
dS = zeros(n,1);
dE = zeros(n,m);
dI = zeros(n,m);
dR = zeros(n,m);

for a = 1:n
    %Calculate the force of infection from each severity level
    FOI = zeros(m,1);
    for s = 1:m
        FOI(s) = para.beta(s)*I(1,s);
    end
   
    dS(a) = -sum(FOI)*S(a);
    
    for s = 1:m
        dE(a,s) = ((1-(1-para.alpha(s))*(1-para.nu(a,s)))*FOI(s)+ para.nu(a,s)*(dot(1-para.alpha,FOI)-(1-para.alpha(s))*FOI(s)))*S(a) - para.epsilon*E(a,s);
        dI(a,s) = para.epsilon*E(a,s) - para.gamma(s)*I(a,s);
        dR(a,s) = para.gamma(s)*I(a,s);
    end    
end

%Convert output to a vector [S; E; I; R]
dPop = zeros((1+3*m)*n,1);
dPop(1:n) = dS;

for ii = 1:m
    dPop(ii*n+1:(ii+1)*n) = dE(:, ii);
    dPop((m+ii)*n+1:(m+ii+1)*n) = dI(:, ii);
    dPop((2*m+ii)*n+1:(2*m+ii+1)*n) = dR(:, ii);
end

end

end
