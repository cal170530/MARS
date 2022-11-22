function w = updateW(numAgents,phi,combMatrix,attacker)

    for i = 1:numAgents
        w(:,i) = [0;0];  %clear weight for each agent
        for j = 1:numAgents
            w(:,i) = w(:,i)+combMatrix(j,i)*phi(:,j);
        end
        if (i == attacker)
            w(:,i) = [2;-2];
        end
    end
                
        
    