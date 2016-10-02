% This MATLAB program checks the feasibility of LMIs from Theorems 1, 2 of the paper 
% A. Selivanov and E. Fridman, "Distributed event-triggered control of diffusion semilinear PDEs," Automatica, vol. 68, pp. 344–351, 2016.
%% Parameters of a chemical reactor 
A=[0 .01; -.45 -.2]; 
B=[1;1]; 
beta=diag([.011, 1.1]); 
D0=diag([.01, .005]); 
Q=diag([10^(-4), 0]);
C=eye(2);
K=[1 0];
rho=.9*ones(2,1);       % quantizer densities 
BC=3;                   % mixed boundary conditions 
Gamma=diag([6 111]);    % the matrix from (5)

%% Point measurements 
l=10;       % domain size
N=25;       % number of sensors
Delta=l/N;  % subdomain size 
delta=2; 
delta1=.9*delta; 

tauM=.0028; epsilon=.09; 
OmegaVal=LMI_Aut16b_th1(A,B,beta,D0,Q,C,K,rho,BC,Gamma,tauM,Delta,epsilon,delta,delta1); 
if ~isempty(OmegaVal)
    syms x
    alpha=eval(solve(x==delta-delta1*exp(x*tauM))); 
    disp(['Theorem 1 (tauM=' num2str(tauM) ', epsilon=' num2str(epsilon) '): Feasible, alpha=' num2str(alpha)]); 
else
    disp(['Theorem 1 (tauM=' num2str(tauM) ', epsilon=' num2str(epsilon) '): Not feasible']); 
end

%% Averaged measurements 
l=10;       % domain size
N=40;       % number of sensors
Delta=l/N;  % subdomain size 
alpha=.3;   % decay rate

tauM=.2859; epsilon=.57; 
OmegaVal=LMI_Aut16b_th2(A,B,beta,D0,Q,C,K,rho,BC,Gamma,tauM,Delta,epsilon,alpha); 
if ~isempty(OmegaVal)
    disp(['Theorem 2 (tauM=' num2str(tauM) ', epsilon=' num2str(epsilon) '): Feasible']); 
else
    disp(['Theorem 2 (tauM=' num2str(tauM) ', epsilon=' num2str(epsilon) '): Not feasible']); 
end