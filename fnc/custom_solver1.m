function [X, chunks, ok] = custom_solver1(allBS, UE, BS_per_km, DEBUG)    
    ok = true; %simple error signaling
    N = max(size(allBS));
    chunks = zeros(N, 1); %this must be a col vector
    for i = 1:N
        chunks(i) = allBS{i}.get_mem_for_UE(BS_per_km);
        if chunks(i) > UE.max_buffer
            chunks(i) = int64(3 * UE.max_buffer / 4);
        end
    end
    
    S = (double(chunks)/1e9) * UE.vel / UE.requested_rate; % [meters]
     
    b = double(chunks)/1e9;
    X = zeros(N,1);
    SX = S' * X;
%     K_old = -1e12;
%     K = -1e11;
    while SX < 1100 && chunks' * X < UE.requested_file_size 
        best_i = 0;
        best_k = -1e12;
%         best_d = zeros(N, 1);
        for i = setdiff(1:N, find(X))
            Y = X;
            Y(i) = 1;              
            
            d = zeros(N, 1);
            indexies = (find(Y))';
            for j = indexies %find(..) finds all nonzeros elements' indexies
                j_x = allBS{j}.pos(1);
                for t = indexies
                    if t~=j
                        t_x = allBS{t}.pos(1);
                        if d(j) >= abs(t_x - j_x) || d(j) == 0
                            d(j) = abs(t_x - j_x);
                        end
                    end
                end
            end       
            
            lambda = 5e2; %tradeoff weight
%             tmp_k = b' * Y - ( ((Y' * Y) \ d' * Y - 1000/(Y' * Y))^2 + (d2' * Y / 2 - 1000)^2 ) / lambda; ((Y' * Y) \ ((d' * Y)) - 1000/(Y' * Y))^2 + 
            tmp_k = b' * Y - ( (d' * Y / 2 - 1000)^2 )/ lambda; %need to divide d'*y by 2 or else count twice the same thing
            if tmp_k > best_k
                best_k = tmp_k;
                best_i = i;
%                 best_d = d;
            end
        end
        
        if DEBUG && best_i > 0
            fprintf('%d ', allBS{best_i}.ID);
        end
            
        if (best_i ~= 0)
            X(best_i) = 1;
%             K_old = K;
%             K = best_k;
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
        title('CUSTOM 1')
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