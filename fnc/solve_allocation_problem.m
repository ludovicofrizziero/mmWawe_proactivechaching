function [X, chunks] = solve_allocation_problem(allBS, UE, DEBUG)
    chunks = zeros(size(allBS)); %this must be a col vector
    for i = 1:size(allBS)
        chunks(i) = allBS{i}.get_mem_for_UE();
    end
    
    f = chunks;
    
    A = chunks'; %This must be a matrix (right now is  1 x n_BS)
    
    b = double(UE.requested_file_size); %this must be a col vector of [(A dot X) x 1] elements
    
    %position whithin x that must be integer values
    intcon = [1:size(allBS)]';
    
    %x must be an integer vector of zeros and ones
    lb = zeros(size(intcon));
    ub = ones(size(intcon));
    
    if DEBUG
        X = intlinprog(-f, intcon, A, b, [], [], lb, ub); %NOTE: -f is required to solve the equivalent minimization problem   
    else
        opts = optimoptions('intlinprog','Display','off');
        X = intlinprog(-f, intcon, A, b, [], [], lb, ub, opts); %NOTE: -f is required to solve the equivalent minimization problem 
    end
    
             
end