function[A]=disk(X,N,Delta);

%% Computes the adjacency matrix for the disk graph
A=zeros(N,N);
for i=1:N;
  for j=1:N;
    if (i~=j & (norm(X(:,i)-X(:,j)) < Delta) );   % Check for proximity
      A(i,j)=1;    
end; end; end;
