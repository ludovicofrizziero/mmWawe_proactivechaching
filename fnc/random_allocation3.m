function [X, chunks, ok] = random_allocation3(allBS, UE, BS_per_km, DEBUG)
    ok = true; %simple error signaling
    N = max(size(allBS));
    chunks = zeros(N, 1); %this must be a col vector
    for i = 1:N
        chunks(i) = allBS{i}.get_mem_for_UE(BS_per_km);
        if chunks(i) > UE.max_buffer
            chunks(i) = int64(3 * UE.max_buffer / 4);
        end
    end
    
    if DEBUG
        tmp1 = double(allBS{1}.get_mem_for_UE(BS_per_km));
        tmp3 = double(UE.max_buffer);
        fprintf('chunks size:\n\torig:\t\t%2.3f GB \n\tue buff:\t%2.3f GB\n', tmp1/8e9, tmp3/8e9);
    end
        
    S = (double(chunks)/1e9) * UE.vel / UE.requested_rate; % [meters]
        
    X = zeros(N,1);
    iter = 1;
    while S' * X < 1100 && chunks' * X < UE.requested_file_size 
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
    
%     if ok
%         I = X' * X; 
%         for i = (find(X))'
%             chunks(i) = int64( (1100 * UE.requested_rate / (I * UE.m_vel)) * 1e9 );
%         end
%     end
    
    %% for debug, plot BS disposition
    if DEBUG
        S = (double(chunks)/1e9) * UE.vel / UE.requested_rate; % [meters]
        figure;
        hold on;
        grid on;
        title('RAND 3');
        for i = 1:N
            plot(allBS{i}.pos(1), allBS{i}.pos(2), '*')
            l = S(i)/2;
            if X(i) > 0
                plot([allBS{i}.pos(1) - l , allBS{i}.pos(1)+S(i)-l], [allBS{i}.pos(2), allBS{i}.pos(2)] - allBS{i}.ID / 10)
            end
            text(allBS{i}.pos(1), allBS{i}.pos(2)+1, strcat('ID: ', int2str(allBS{i}.ID)));
%             text(allBS{i}.pos(1), allBS{i}.pos(2)+2, strcat('w: ', num2str(n(i)*w(i)))); %weigth of BS
            if X(i) > 0
                text(allBS{i}.pos(1), allBS{i}.pos(2)+3, 'yes') %num2str(double(allBS{i}.memory), '%1.3e'));
            end
        end
        hold off;
        
        fprintf('SX: %f\n', S' * X);
    end
    %%
end