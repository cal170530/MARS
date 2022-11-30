%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear;
clc;

figure
%% Proximity disk
Delta=3; % interaction distance

%% You should try your solutions on three different, initial networks that
%% you load using the following function call


% n = dimension of each robot (n>1)
% N = total number of robots
% X = n x N vector containing the initial robot positions

%mode = 'static';
%mode = 1 %mode = 'moving';
mode = 1;
%% Some numerical integration parameters
dt=0.003; % numerical steplength
Tf=4; % final time
network_no=4;


if mode == 2
   network_no=2; 
   [X,n,N]=load_network(network_no,Delta);
   numAgents = N-2; 
   w0 = X(:,N);
   w1 = X(:,N-1);
   attacker = N-2;
   groupDivision = 5;
end
if mode == 1
   network_no=4; 
   [X,n,N]=load_network(network_no,Delta);
   numAgents = N-1; 
   attacker = N-1;
   w0 = X(:,N); %Estimation Parameter: Nth X position used as a static parameter
   groupDivision = 100;
end
t=0; 
iter=2;
MSD = 0;

dist = 0;
p = 3;
%%%%%%Adaptive Weight Variables%%%%%%%%%%
combMatrix = zeros(numAgents);  %%nxn matrix of adaptive weights a(j,i), left stochastic
phi = zeros(2,1,numAgents); %% 2x1 vectors for each agent's intermediate parameter estimation 
gamma = zeros(numAgents);   %% nxn matrix of variance estimates of agent j according to agent i, gamma(j,i)
%Initialize regression sequence and weights for agents
for k = 1:numAgents
u(:,k) = mvnrnd(0,.8,2000);
w(:,k)= [normrnd(0,1);normrnd(0,1)];
vel(:,k) = [0;0];
groupVel(:,k) = [0;0];
del(:,k) = [0;0];
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
step = 0.05;
dist = @(xi, xj) sqrt( (xi - xj)' * (xi - xj));
%weight = @(xi, xj) weightfcn(dist(xi, xj)); 
% H for hysterisys
epsilon = .05;
H = zeros(size(disk(X,N,Delta)));
dmin = inf;


while (t<=Tf);

%% A is the adjacency matrix associated with the system
%% using a disk graph. 
A=disk(X,N,Delta);
DX=zeros(n,N); %% Here is where we store the derivatives 
w0 = X(:,N);
for i=1:numAgents
    for j=1:numAgents;
      if (A(i,j)==1);
        dmin = min(dist(X(:,j),X(:,i)), dmin);
        if H(i, j) == 0     % this is a new edge
             if dist(X(:,j),X(:,i)) < Delta-epsilon  % don't head towards
                H(i, j) = 1;
             end
        end
      end
    end
%%%%%Adapt Step%%%%%%
if (i ~= attacker)
ux = u(iter:-1:iter-1,i)'; % Current regression vector
v = normrnd(0,.15);     %Gaussian noise
if i<groupDivision
d = ux*w0+v;
else
d = ux*w1+v;
end
phi(:,i) = w(:,i) +step*ux'*(d - ux*w(:,i)); 
end
if (i == attacker)
phi(:,i) = [2;-2];
end
%%%%%%%%%%%%%%%%%%%%%  
end;
combMatrix = updateCombinationMatrix(combMatrix,numAgents,w,phi,gamma,A,attacker);
w = updateW(numAgents,phi,combMatrix,attacker);
groupVel = ComputeGroupVelocity(combMatrix,w,vel,groupVel,numAgents);
del = ComputeDel(A,numAgents,X);

%% Update the states using an Euler approximation
  for i=1:numAgents
   % DX(:,i)=DX(:,i)+(w(:,i)-X(:,i));
   if (i ~= attacker)
   vel(:,i) = .5*(w(:,i)-X(:,i)) + .5*groupVel(:,i)+del(:,i) 
   X(:,i)=X(:,i)+dt.*vel(:,i);
   end
  end
if mode == 1
X(:,N) = [-2;9] +2*[cos(t);-sin(t)];
end
%% Update time
  t=t+dt;


%% Plot the solution every 10 iterations
  if (mod(iter,10)==0);
    plotsol(X,N,A,Delta,attacker);
  end;

  iter=iter+1;
end;

for k = 1:numAgents
    if k~= attacker
        MSD = MSD + norm(w0-w(:,k))^2;
    end
end
MSD = MSD/numAgents;
MSD


