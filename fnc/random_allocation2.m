function [X, chunks, ok] = random_allocation2(allBS, UE, BS_per_km, DEBUG)
    ok = true; %simple error signaling
    N = max(size(allBS));
    chunks = zeros(N, 1); %this must be a col vector
    for i = 1:N
        chunks(i) = int64(UE.max_buffer * rand(1)); 
    end
        
    S = (double(chunks)/1e9) * UE.vel / UE.requested_rate; % [meters]
    
    for i = 1:N
        if chunks(i) > UE.max_buffer
            chunks(i) = 0;
        end
    end
        
    X = zeros(N,1);
    iter = 1;
    while S' * X < 1000
        i = randi(N);
        X(i) = 1;   
        iter = iter + 1;
        if iter == N
            break; %to avoid bugs
        end
    end      
    
    if S' * X < 1000
        ok = false;
    end   
    
    %% for debug, plot BS disposition
    if DEBUG
        figure;
        hold on;
        title('RAND 2')
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