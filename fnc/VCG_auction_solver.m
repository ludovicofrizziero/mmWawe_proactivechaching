function [X, chunks, ok] = VCG_auction_solver(allBS, UE, BS_per_km, DEBUG)
%     ok = true; %simple error signaling
%     N = max(size(allBS));
%     chunks = zeros(N, 1); %this must be a col vector
%     for i = 1:N
%         chunks(i) = allBS{i}.get_mem_for_UE(BS_per_km);
%     end
%         
%     S = (double(chunks)/1e9) * UE.vel / UE.requested_rate; % [meters]
%     
%     w = 1/N * (S - mean(S))*(S - mean(S))'; % cov matrix
%     %w = 1/N * (S - 1000/BS_per_km)*(S - 1000/BS_per_km)'; % 'look alike' cov matrix
%     w = diag(w);
%     w = max(1e-9, w);
%     
%     for i = 1:N
%         if chunks(i) > UE.max_buffer
%             chunks(i) = 0;
%         end
%     end
%     
%     K = 0;
%     K_old = -1;
%     I = 1:N;
%     X = [];
%     i = 0;
%     while K > K_old && ~isempty(I)
%         [~, i] = max(chunks(I)); %argmax
%         X = [X; I(i)];
%         I(i) = [];
%         
%        
%         K_old = K; 
%         B = sum(double(chunks(X)) / 1e9);
%         SX = sum(S(X));        
%         a = ((double(chunks(X)) ./ 1e9) ./ B);
%         b = (1 ./ w(X));
%         c = exp(-(SX-1000)^2 / (2 * 1e3));
%         K = sum( a .* b ) * c;        
%     end
% 
%     X = setdiff(X, i); %remove last entry since caused a diminishing of K
%     
%     tmp = zeros(N,1);
%     tmp(X) = 1;
%     X = tmp;
%     
%     
%     %% for debug, plot BS disposition
%     if DEBUG
%         S = (double(chunks)/1e9) * UE.vel / UE.requested_rate;
%         figure;
%         hold on;
%         grid on;
%         for i = 1:N
%             plot(allBS{i}.pos(1), allBS{i}.pos(2), '*')
%             l = S(i)/2;
%             plot([allBS{i}.pos(1) - l , allBS{i}.pos(1)+S(i)-l], [allBS{i}.pos(2), allBS{i}.pos(2)] - allBS{i}.ID / 10)
%             text(allBS{i}.pos(1), allBS{i}.pos(2)+1, strcat('ID: ', int2str(allBS{i}.ID)));
%             text(allBS{i}.pos(1), allBS{i}.pos(2)+2, strcat('w: ', num2str(w(i)))); %weigth of BS
%             if X(i) > 0
%                 text(allBS{i}.pos(1), allBS{i}.pos(2)+3, 'yes')
%             end
%         end
%         hold off;
%         grid off;
%         fprintf('SX: %f\n', S' * X);
%     end
%     %%


    ok = true; %simple error signaling
    N = max(size(allBS));
    chunks = zeros(N, 1); %this must be a col vector
    for i = 1:N
        chunks(i) = int64(rand(1) * UE.max_buffer); %allBS{i}.get_mem_for_UE(BS_per_km);
        if chunks(i) > UE.max_buffer
            chunks(i) = 0;
        end
    end
    
    S = (double(chunks)/1e9) * UE.vel / UE.requested_rate; % [meters]
    
    c = 1000/BS_per_km; 
    b = double(chunks)/1e9;
    X = zeros(N,1);
    SX = S' * X;
    K_old = -1e12;
    K = -1e11;
    while SX < 1100 && chunks' * X < UE.requested_file_size 
        best_i = 0;
        best_k = -1e12;

        for i = setdiff(1:N, find(X))
            Y = X;
            Y(i) = 1;
                        
            d2 = zeros(N, 1);
            for j = (find(Y))' %find(..) finds all nonzeros elements' indexies
                j_x = allBS{j}.pos(1);
                for t = (find(Y))'
                    if t~=j
                        t_x = allBS{t}.pos(1);
                        if d2(j) >= abs(t_x - j_x) || d2(j) == 0
                            d2(j) = abs(t_x - j_x);
                        end
                    end
                end
            end       
            
            lambda = 5e4; %tradeoff weight
%             tmp_k = b' * Y - ( ((Y' * Y) \ d' * Y - 1000/(Y' * Y))^2 + (d2' * Y / 2 - 1000)^2 ) / lambda; ((Y' * Y) \ ((d' * Y)) - 1000/(Y' * Y))^2 + 
            tmp_k = b' * Y - ( (d2' * Y / 2 - 1000)^2 )/ lambda; %need to divide d'*y by 2 or else count twice the same thing
            if tmp_k > best_k
                best_k = tmp_k;
                best_i = i;
%                 best_d = d;
            end
        end
        
        if DEBUG
            fprintf('%d ', allBS{best_i}.ID);
        end
            
        if (best_i ~= 0)
            X(best_i) = 1;
            K_old = K;
            K = best_k;
            SX = S' * X;
        else
            %disp('problem -> maximum reached before constraints were respected.')
            ok = false;
            break;
        end
    end

    
    %% for debug, plot BS disposition
    if DEBUG
        figure;
        hold on;
        grid on;
        for i = 1:N
            plot(allBS{i}.pos(1), allBS{i}.pos(2), '*')
            l = S(i)/2;
            plot([allBS{i}.pos(1) - l , allBS{i}.pos(1)+S(i)-l], [allBS{i}.pos(2), allBS{i}.pos(2)] - allBS{i}.ID / 10)
            text(allBS{i}.pos(1), allBS{i}.pos(2)+1, strcat('ID: ', int2str(allBS{i}.ID)));
%             text(allBS{i}.pos(1), allBS{i}.pos(2)+2, strcat('w: ', num2str(n(i)*w(i)))); %weigth of BS
            if X(i) > 0
                text(allBS{i}.pos(1), allBS{i}.pos(2)+3, 'yes') %num2str(double(allBS{i}.memory), '%1.3e'));
            end
        end
        hold off;
        grid off;
        
        fprintf('SX: %f\n', S' * X);
    end
    %%



end