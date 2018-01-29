function [X, chunks] = VCG_auction_solver(allBS, UE, DEBUG)

    N = max(size(allBS));
    chunks = zeros(N, 1); %this must be a col vector
    for i = 1:N
        chunks(i) = allBS{i}.get_mem_for_UE();
    end
        
    S = (double(chunks)/1e9) * UE.vel / UE.requested_rate; % [meters]
    w = ones(N, 1);
    for i = 1:N
        l = abs(allBS{i}.pos(2) - allBS{i}.sharedData.UE.pos(2)) * sin(pi/3) / sin(pi/6);
        for j = 1:N %BS are not ordered by distance hence need to serch all of them each time
            if  i ~= j
                if allBS{j}.pos(1) >= allBS{i}.pos(1) - l && allBS{j}.pos(1) <= allBS{i}.pos(1) + S(i) - l
                    w(j) = w(j) + 1.0;
                end
            end
        end
    end
    
    for i = 1:N
        if chunks(i) > UE.max_buffer
            chunks(i) = 0;
        end
    end
    
    K = 0;
    K_old = -1;
    I = 1:N;
    X = [];
    i = 0;
    %perm = randperm(N);
    %chunks  = chunks(perm); %this is needed only to avoid selecting in order only the first BSs when they all have same chunks
    while K > K_old && ~isempty(I)
        [~, i] = max(chunks(I)); %argmax
        X = [X; I(i)];
        I(i) = [];
        
       
        K_old = K; 
        B = sum(double(chunks(X)) / 1e9);
        SX = sum(S(X));        
        a = ((double(chunks(X)) ./ 1e9) ./ B);
        b = (1 ./ w(X));
        c = exp(-(SX-1000)^2 / (2 * 1e3));
        K = sum( a .* b ) * c;        
    end

    X = setdiff(X, i); %remove last entry since caused a diminishing of K
    
    tmp = zeros(N,1);
    tmp(X) = 1;
    X = tmp;
    
    
    %% for debug, plot BS disposition
    if DEBUG
        figure;
        hold on;
        for i = 1:N
            plot(allBS{i}.pos(1), allBS{i}.pos(2), '*')
            l = abs(allBS{i}.pos(2) - allBS{i}.sharedData.UE.pos(2)) * sin(pi/3) / sin(pi/6);
            plot([allBS{i}.pos(1) - l , allBS{i}.pos(1)+S(i)-l], [allBS{i}.pos(2), allBS{i}.pos(2)] - allBS{i}.ID / 10)
            text(allBS{i}.pos(1), allBS{i}.pos(2)+1, strcat('ID: ', int2str(allBS{i}.ID)));
            text(allBS{i}.pos(1), allBS{i}.pos(2)+2, strcat('w: ', num2str(w(i)))); %weigth of BS
            if X(i) > 0
                text(allBS{i}.pos(1), allBS{i}.pos(2)+3, 'yes') %num2str(double(allBS{i}.memory), '%1.3e'));
            end
        end
        hold off;
        
        fprintf('SX: %f\n', S' * X);
    end
    %%
end