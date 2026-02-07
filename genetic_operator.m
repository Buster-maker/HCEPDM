function f  = genetic_operator(parent_chromosome, M, V, mum, l_limit, u_limit, item, tao,cntao,probID)

%% function f  = genetic_operator(parent_chromosome, M, V, mu, mum, l_limit, u_limit)
% 
% This function is utilized to produce offsprings from parent chromosomes.
% The genetic operators corssover and mutation which are carried out with
% slight modifications from the original design. For more information read
% the document enclosed. 
%
% parent_chromosome - the set of selected chromosomes.
% M - number of objective functions
% V - number of decision varaiables
% mu - distribution index for crossover (read the enlcosed pdf file)
% mum - distribution index for mutation (read the enclosed pdf file)
% l_limit - a vector of lower limit for the corresponding decsion variables
% u_limit - a vector of upper limit for the corresponding decsion variables
%
% The genetic operation is performed only on the decision variables, that
% is the first V elements in the chromosome vector. 




[N,m] = size(parent_chromosome);

clear m


interf=[];
p = 1;
for i = 1 :3: N
     if rand(1) < 0.9 
     % Initialize the children to be null vector.
     child_1 = [];
             
     % Get the chromosome information for each randomnly selected
     % parents
     parent_1 = parent_chromosome(i,1:V);
     parent_2 = parent_chromosome(i+1,1:V);
     parent_3 = parent_chromosome(i+2,1:V);
     % Perform corssover for each decision variable in the chromosome.
       for j = 1 : V
            
         
            % SBX (Simulated Binary Crossover).
            % For more information about SBX refer the enclosed pdf file.
            % Generate a random number
                    
            
             child_1(j) =parent_1(j) + 0.5*(parent_2(j)-parent_3(j));
             
            
            
            
            % Make sure that the generated element is within the specified
            % decision space else set it to the appropriate extrema.
            if child_1(j) > u_limit(j)
                child_1(j) = u_limit(j);
            elseif child_1(j) < l_limit(j)
                child_1(j) = l_limit(j);
            end
            
       end
%        
%                apha=1;
%         I=abs(parent_1-parent_2);
%         L1=max(l_limit,parent_1-I*apha);
%         U1=min(u_limit,parent_1+I*apha);
%         L2=max(l_limit, parent_2-I*apha);
%         U2=min(u_limit, parent_2+I*apha);
%         r1=rand(1);
%         child_1=L1+r1*(U1-L1);           
%         r2=rand(1);
%         child_2=L2+r2*(U2-L2); 
       
       interf(p,:)=child_1;
       p=p+1;
     else
        interf(p,:) = parent_chromosome(i,1:V);
        p=p+1;
     end
end

        
      
for i=1:p-1
    if rand(1)<0.1
        
        child_3 = interf(i,:);
        % Perform mutation on eact element of the selected parent.
        for j = 1 : V
           r(j) = rand(1);
           if r(j) < 0.5
               delta(j) = (2*r(j))^(1/(mum+1)) - 1;
           else
               delta(j) = 1 - (2*(1 - r(j)))^(1/(mum+1));
           end
           % Generate the corresponding child element.
           child_3(j) = child_3(j) + delta(j);
           % Make sure that the generated element is within the decision
           % space.
           if child_3(j) > u_limit(j)
               child_3(j) = u_limit(j);
           elseif child_3(j) < l_limit(j)
               child_3(j) = l_limit(j);
           end
        end
        interf(i,:)=child_3;
    end
end
p;
for i=1:p-1
    
    interf(i,V + 1: M + V) = evaluate_objective(interf(i,:), M, V,item,tao,cntao,probID);
end
f=interf;
