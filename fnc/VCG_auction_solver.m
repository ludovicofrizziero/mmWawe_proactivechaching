function [X, chunks] = VCG_auction_solver(allBS, UE, BS_per_km, DEBUG)

    N = max(size(allBS));
    chunks = zeros(N, 1); %this must be a col vector
    for i = 1:N
        chunks(i) = allBS{i}.get_mem_for_UE(BS_per_km);
    end
        
    S = (double(chunks)/1e9) * UE.vel / UE.requested_rate; % [meters]
    
    w = 1/N * (S - mean(S))*(S - mean(S))'; % cov matrix
    %w = 1/N * (S - 1000/BS_per_km)*(S - 1000/BS_per_km)'; % 'look alike' cov matrix
    w = diag(w);
    w = max(1e-9, w);
    
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
        S = (double(chunks)/1e9) * UE.vel / UE.requested_rate;
        figure;
        hold on;
        grid on;
        for i = 1:N
            plot(allBS{i}.pos(1), allBS{i}.pos(2), '*')
            l = S(i)/2;
            plot([allBS{i}.pos(1) - l , allBS{i}.pos(1)+S(i)-l], [allBS{i}.pos(2), allBS{i}.pos(2)] - allBS{i}.ID / 10)
            text(allBS{i}.pos(1), allBS{i}.pos(2)+1, strcat('ID: ', int2str(allBS{i}.ID)));
            text(allBS{i}.pos(1), allBS{i}.pos(2)+2, strcat('w: ', num2str(w(i)))); %weigth of BS
            if X(i) > 0
                text(allBS{i}.pos(1), allBS{i}.pos(2)+3, 'yes')
            end
        end
        hold off;
        grid off;
        fprintf('SX: %f\n', S' * X);
    end
    %%
end