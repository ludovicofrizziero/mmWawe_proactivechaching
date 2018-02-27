function [X, chunks] = solve_allocation_problem(allBS, UE, DEBUG)   
    N = max(size(allBS));
    chunks = zeros(N, 1); %this must be a col vector
    for i = 1:N
        chunks(i) = allBS{i}.get_mem_for_UE();
    end
        
    S = (double(chunks)/1e9) * UE.vel / UE.requested_rate; % [meters]
    f = zeros(N, 1);
    %f = ones(N, 1);
    %g = ones(N, 1);
    for i = 1:N
        l = abs(allBS{i}.pos(2) - allBS{i}.sharedData.UE.pos(2)) * sin(pi/3) / sin(pi/6);
        for j = 1:N %BS are not ordered by distance hence need to serch all of them each time
            if  i ~= j
                if allBS{j}.pos(1) >= allBS{i}.pos(1) - l && allBS{j}.pos(1) <= allBS{i}.pos(1) + S(i) - l
                    f(j) = f(j) + 1.0;
                    %f(i) = f(i) + 1.0;
                end
            end
        end
    end
    
%     for i = 1:N
%         for j = 1:N %BS are not ordered by distance hence need to serch all of them each time
%             if  i ~= j
%                 if allBS{j}.pos(1) >= allBS{i}.pos(1) - S(i)/4 && allBS{j}.pos(1) <= allBS{i}.pos(1) + 3*S(i)/4
%                     g(i) = g(i) + f(j);
%                 end
%             end
%         end
%         g(i) = g(i) + f(i);
%     end
%     
%     f = g;
    
    % A must be a [K x len(X)] matrix
    A = chunks'; 
    A = [A; diag(chunks)];
    A = [A; - S'];
    
    % b must be a [K x 1] vector (that is, such that A*X <= b is determined)
    b = double(UE.requested_file_size); 
    b = [b; ones(N, 1) * double(UE.max_buffer)];
    b = [b; -1000];
    
    %position whithin x that must be integer values
    intcon = (1:size(allBS))';
    
    %x must be an integer vector of zeros and ones
    lb = zeros(size(intcon));
    ub = ones(size(intcon));
    
    if DEBUG
        X = intlinprog(f, intcon, A, b, [], [], lb, ub); 
    else
        opts = optimoptions('intlinprog','Display','off');
        X = intlinprog(f, intcon, A, b, [], [], lb, ub, opts); 
    end
    
    %% for debug, plot BS disposition
    if DEBUG
        figure;
        hold on;
        for i = 1:N
            plot(allBS{i}.pos(1), allBS{i}.pos(2), '*')
            l = abs(allBS{i}.pos(2) - allBS{i}.sharedData.UE.pos(2)) * sin(pi/3) / sin(pi/6);
            plot([allBS{i}.pos(1) - l , allBS{i}.pos(1)+S(i)-l], [allBS{i}.pos(2), allBS{i}.pos(2)])
            text(allBS{i}.pos(1), allBS{i}.pos(2)+1, strcat('ID: ', int2str(allBS{i}.ID)));
            text(allBS{i}.pos(1), allBS{i}.pos(2)+2, strcat('w: ', num2str(f(i)))); %weigth of BS
            if X(i) > 0
                text(allBS{i}.pos(1), allBS{i}.pos(2)+3, 'yes') %num2str(double(allBS{i}.memory), '%1.3e'));
            end
        end
        hold off;
    end
    %%
    
             
end