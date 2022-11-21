function w = UpdateWeight(numAgents,gamma,A)
    for i = 1:numAgents
        numNeighbors = sum(A(i,:))+1;
        w(:,i) = gamma(:,i)/numNeighbors;
        for j = 1:numAgents
            if A(i,j) == 1
                w(:,i) = w(:,i) + gamma(:,j)/numNeighbors;
            end
        end
    end
                
        
    