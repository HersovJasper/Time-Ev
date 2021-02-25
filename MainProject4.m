%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Investigates the arrow of time

clear all
tic

global H
hbar = 1;

%%%%%%%%%%%%%%%%%%
%% Parameters	%%
%%%%%%%%%%%%%%%%%%
N=3;        % Number of atoms
M=5;        % Number of wells/modes - fixed to 2
E = zeros(M,1);   % Energy of well set to zero
%E=rand(M,1);      % could look at random well depths to see if this has a difference 
J  = 1;    % Tunnelling constant
U  = 0.1;    % Interaction constant (-ve means attractive)

%*********** Load Matrices ***********************
% Loads matrices and if they haven't been created it makes them first
data=sprintf('BHM_data_N=%i_M=%i/Hamiltonian',N,M);
dir=sprintf('BHM_data_N=%i_M=%i',N,M);

if exist(dir,'dir')
  load(data);
else
  ham = func_ham_maker_BHM(N,M);
  load(data);
end
L=length(basis);

%%%%%%%%%%%%%% Hamiltonian %%%%%%%%%%%%%%%%%%%%		
HET = zeros(L,L);	
for m=1:M
  HET = HET + E(m)*HE(:,:,m);
endfor

if N>1
  H = J*HJ + U*HU + HET;
else
  H = J*HJ + HET;
end

[vec, val] = eig(H);



%%%%%%%%%%% Time Evolution %%%%%%%%%%%%%%%%%%%%%%%%
%********** Initial state ************************
q=1;
for m=1:L
  for p=1:M
    if basis(m,p)==N
      noon(q)=m;
      q=q+1;
    end
  end
end

psi = zeros(L,1);
%psi(noon(round(M/2)))=1;
psi(1)=1;

%*********** Evolution Parameters *****************
tmax = 20;
dt = 0.1;
t = 0:dt:tmax;
Lt = length(t);
psiE = vec'*psi;


Exp = zeros(Lt,M);
P_noon = zeros(Lt,M);
P_all = zeros(Lt,L);
       


for n=1:Lt
  psiEt = exp(-i*diag(val)*t(n)/hbar).*psiE;
  psit = vec*psiEt;
  
  for m=1:M
    Exp(n,m)=psit'*diag(basis(:,m))*psit;
    P_noon(n,m) = abs(psit(noon(m)))^2;
    P_all(n,:) = abs(psit).^2;
  end
  
end

figure(1)
plot(t,Exp)

figure(2)
plot(t,P_noon)

figure(3)
plot(t,P_all)

figure(4)
image(P_all*100)

figure(5)
image(Exp*50)


toc
