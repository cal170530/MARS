%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear;
clc;

figure
%% Proximity disk
Delta=3; % interaction distance

%% You should try your solutions on three different, initial networks that
%% you load using the following function call

network_no=4; %% or 2 or 3
[X,n,N]=load_network(network_no,Delta);
% n = dimension of each robot (n>1)
% N = total number of robots
% X = n x N vector containing the initial robot positions


%% Some numerical integration parameters
dt=0.003; % numerical steplength
Tf=3; % final time
numAgents = N-1;
t=0; 
iter=2;
MSD = 0;
dist = 0;
p = 3;
gamma = zeros(2,1,numAgents);

%weightfcn = @(d) tan(-pi/2 + pi * d/Delta) + b;
 weightfcn = @(d) (.5 .* (2 .* Delta - d)./(Delta - d).^2) * (d - dist);
dist = @(xi, xj) sqrt( (xi - xj)' * (xi - xj));
weight = @(xi, xj) weightfcn(dist(xi, xj)); 

% H for hysterisys
epsilon = .05;
H = zeros(size(disk(X,N,Delta)));
dmin = inf;
w = [];
%Initialize regression sequence and weights for agents
for k = 1:numAgents
u(:,k) = mvnrnd(0,.8,2000);
w(:,k)= [normrnd(0,1);normrnd(0,1)];
end
%Estimation Parameter
w0 = X(:,N);

step = 0.05;

while (t<=Tf);

%% A is the adjacency matrix associated with the system
%% using a disk graph. 
A=disk(X,N,Delta);
DX=zeros(n,N); %% Here is where we store the derivatives 
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
ux = u(iter:-1:iter-1,i)';
d = ux*w0+normrnd(0,.15);
gamma(:,i) = gamma(:,i) +step*ux'*(d - ux*w(:,i)); 
%%%%%%%%%%%%%%%%%%%%%  
end;  
w = UpdateWeight(numAgents,gamma,A);
%% Update the states using an Euler approximation
  for i=1:numAgents
    DX(:,i)=DX(:,i)+(w(:,i)-X(:,i));
    X(:,i)=X(:,i)+dt.*DX(:,i);
  end

%% Update time
  t=t+dt;


%% Plot the solution every 10 iterations
  if (mod(iter,10)==0);
    plotsol(X,N,A,Delta);
  end;

  iter=iter+1;
end;

for k = 1:numAgents
    MSD = MSD + norm(w0-w(:,k))^2;
end
MSD = MSD/numAgents;
MSD


