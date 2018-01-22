function solve_allocation_problem(allBS, UE)
    tmp = zeros(size(allBS, 1)); %this must be a col vector
    for i = 1:size(allBS)
        tmp(i) = allBS{i}.get_mem_for_UE();
    end
    
    f = tmp;
    
    A = tmp'; %This must be a matrix (right now is  1 x n_BS)
    
    b = UE.requested_file_size; %this must be a col vector of [(A dot X) x 1] elements
    
    %position whithin x that must be integer values
    intcon = 1:size(allBS);
    
    %x must be an integer vector of zeros and ones
    lb = zeros(size(intcon), 1);
    ub = ones(size(intcon), 1);
    
    X = intlinprog(-f, intcon, A, b, [], [], lb, ub); %NOTE: -f is required to solve the equivalent minimization problem    
    
    for i = 1:size(X)
        allBS{i}.allocate_memory_for_ue(tmp(i) * X(i));
    end
end