function del = ComputeDel(A,numAgents,X)
dist = @(xi, xj) sqrt( (xi - xj)' * (xi - xj));
tau = 1.5;
for k = 1:numAgents
del(:,k) = [0;0];
neighbors = sum(A(k,:));
b = 1/(neighbors+1)
for j = 1:numAgents
    if A(k,j) == 1
        del(:,k) = del(:,k)+ b*(dist(X(:,j),X(:,k))- tau)*(X(:,j)-X(:,k))/dist(X(:,j),X(:,k))
    end
end
end
