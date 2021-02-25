%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Main.m
%-------------------------------------------------------------------------------
% Loads Hamiltonian matrices for chosen N and M values

clear all
tic

global H
hbar = 1;

%%%%%%%%%%%%%%%%%%
%% Parameters	%%
%%%%%%%%%%%%%%%%%%
N=5;M=4;
data=sprintf('data_N=%i_M=%i/Hamiltonian',N,M);
load(data)
L=length(basis);S=size(basis);
g=0.d0; % Coupling between environment and system

%************** Thermal atoms ******************%
EL0 = 0; ER0 =0 ; %	Energy of lower levels
EL1 = 1; ER1 =1; % Energy of upper levels
J0  = 0.15; J1  = 0.2; % Tunneling between wells
U0  = 2/N; UL0=U0; UR0=U0; % Interaction strength in 0 level
U1  = 1/N; UL1=U1; UR1=U1; % Interaction strength in 0 level
UL01 = 0.1; UR01 = 0.1; % Level switching interaction
%***********************************************%

%************** Single atom ********************%
T=0.1; % Tunneling of the single atom
gL0 = g; gR0 = g; % Interaction between single atom and ground levels
gL1 = g; gR1 = g; % Interaction between single atom and excited levels
%***********************************************%


%************ Initial Wavefunction *************%
% Thermal state so 80% in ground levels and 20% in excited levels
% Put more atoms in one well
N0 = floor(0.8*N); N1 = N - floor(0.8*N);
conf=[floor(0.9*N0),N0-floor(0.9*N0),floor(0.9*N1),N1-floor(0.9*N1)];

% Initial configurations of the system
conf0=[conf,1,0];							
conf1=[conf,0,1];

% Find states with these configurations
ind0=find(sum(abs(basis-ones(L,1)*conf0),2)==0);
ind1=find(sum(abs(basis-ones(L,1)*conf1),2)==0);

% Initial wavefunction
psi=zeros(L,1); psi(ind0)=1/sqrt(2); psi(ind1)=1/sqrt(2);

%************ Hamiltonian **********************	
% Individual saved Hamiltonian matrices:
% 'HU0L','HU0R','HU1L','HU1R','HU01L','HU01R',
% 'HE0L','HE0R','HE1L','HE1R',
% 'HJ0','HJ1',
% 'Hg0L','Hg0R','Hg1L','Hg1R',
% 'HT'		
H = EL0*HE0L + ER0*HE0R + EL1*HE1L + ER1*HE1R...
    + J0*HJ0 + J1*HJ1...
    + UL0*HU0L + UR0*HU0R + UL1*HU1L + UR1*HU1R...
    + UL01*HU01L + UR01*HU01R...
    + T*HT...
    + gL0*Hg0L + gR0*Hg0R + gL1*Hg1L + gR1*Hg1R;

[vec, val] = eig(H);
  

%%%%%%%%%%%%%%%%%%
%% Evolution	%%
%%%%%%%%%%%%%%%%%%  
%*********** Evolution Parameters *****************
tmax = 10;
dt = 1;
t = 0:dt:tmax;
Lt = length(t);
psiE = vec'*psi;

nL0Exp = zeros(Lt,1);
nR0Exp = zeros(Lt,1);
nL1Exp = zeros(Lt,1);
nR1Exp = zeros(Lt,1);


for n=1:Lt
  psiEt = exp(-i*diag(val)*t(n)/hbar).*psiE;
  psit = vec*psiEt;
  nL0Exp(n)=psit'*diag(basis(:,1))*psit;
  nR0Exp(n)=psit'*diag(basis(:,2))*psit;
  nL1Exp(n)=psit'*diag(basis(:,3))*psit;
  nR1Exp(n)=psit'*diag(basis(:,4))*psit;
endfor

plot(t,nL0Exp,t,nR0Exp,t,nL1Exp,t,nR1Exp,t,nL0Exp+nR0Exp+nL1Exp+nR1Exp)






toc
