function [X, chunks] = random_allocation(allBS, UE, BS_per_km, DEBUG)

    N = max(size(allBS));
    chunks = zeros(N, 1); %this must be a col vector
    for i = 1:N
        chunks(i) = allBS{i}.get_mem_for_UE(BS_per_km);
    end
        
    S = (double(chunks)/1e9) * UE.vel / UE.requested_rate; % [meters]
    
    for i = 1:N
        if chunks(i) > UE.max_buffer
            chunks(i) = 0;
        end
    end
        
    SX = 0;
    I = 1:N;
    X = [];
    while SX < 1100 && ~isempty(I)
        i = randi(size(I));
        X = [X; I(i)];
        I(i) = [];
        SX = sum(S(X));
    end
   
    tmp = zeros(N,1);
    tmp(X) = 1;
    X = tmp;
    
    
    %% for debug, plot BS disposition
    if DEBUG
        figure;
        hold on;
        for i = 1:N
            plot(allBS{i}.pos(1), allBS{i}.pos(2), '*')
            l = S(i)/2;
            plot([allBS{i}.pos(1) - l , allBS{i}.pos(1)+S(i)-l], [allBS{i}.pos(2), allBS{i}.pos(2)] - allBS{i}.ID / 10)
            text(allBS{i}.pos(1), allBS{i}.pos(2)+1, strcat('ID: ', int2str(allBS{i}.ID)));
            %text(allBS{i}.pos(1), allBS{i}.pos(2)+2, strcat('w: ', num2str(w(i)))); %weigth of BS
            if X(i) > 0
                text(allBS{i}.pos(1), allBS{i}.pos(2)+3, 'yes') %num2str(double(allBS{i}.memory), '%1.3e'));
            end
        end
        hold off;
        
        fprintf('SX: %f\n', S' * X);
    end
    %%
end