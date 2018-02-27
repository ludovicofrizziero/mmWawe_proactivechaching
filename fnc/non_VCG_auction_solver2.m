function [X, chunks, ok] = non_VCG_auction_solver2(allBS, UE, BS_per_km, DEBUG)
    ok = true; %simple error signaling
    N = max(size(allBS));
    chunks = zeros(N, 1); %this must be a col vector
    for i = 1:N
        chunks(i) = allBS{i}.get_mem_for_UE(BS_per_km);
    end
            
    d = zeros(N,1);
    for i = 1:N
        bsi = allBS{i};
        tmp = zeros(N, 1);
        for j = 1:N
            if i ~= j
                bsj = allBS{j};
                tmp(j) = bsi.pos(1) - bsj.pos(1);                
            end
        end
        tmp = abs(tmp);
        tmp = sort(tmp);
        d(i) = mean(tmp(1:4));
    end
    
    S = (double(chunks)/1e9) * UE.vel / UE.requested_rate; % [meters]
    w1 = 1/N * (S - mean(S))*(S - mean(S))'; % cov matrix
    w1 = diag(w1);
    w1 = max(1e-9, w1);
    w2 = 1/N * (d - 1000/BS_per_km)*(d - 1000/BS_per_km)'; % 'look alike' cov matrix
    w2 = diag(w2);
    w2 = max(1e-9, w2);
    w = w1 .* w2;
    
    for i = 1:N
        if chunks(i) > UE.max_buffer
            chunks(i) = 0;
        end
    end
    
    K = 0;
    K_old = -1;
    I = 1:N;
    X = []; 
    n = [];
    while K > K_old && ~isempty(I)    
        n = ones(N,1);
%         for i = I(X)
%             l1 = S(i)/2;
%             start1 = allBS{i}.pos(1) - l1;
%             stop1 = allBS{i}.pos(1) + S(i) - l1;
%             for j = setdiff(I, X) 
%                 l2 = S(j)/2;
%                 start2 = allBS{j}.pos(1) - l2;
%                 stop2 = allBS{j}.pos(1) + S(j) - l2;
%                 
%                 if start2 >= start1 && start2 <= stop1 || stop2 <= stop1 && stop2 >= start1
%                     n(j) = n(j) + 1;
%                 end
% %                 if allBS{j}.pos(1) < 1000/BS_per_km || allBS{j}.pos(1) > 1000 - 1000/BS_per_km
% %                     n(j) = n(j) + 1; % otherwise BS on the border are advantaged
% %                 end
% %                 if  abs(allBS{j}.pos(1) -  allBS{i}.pos(1)) < (1000/ BS_per_km) / 2
% %                     n(j) = n(j) + 1;
% %                 end
%             end
%         end
        n = ones(N,1);
        Z_old = 0;
        best = 0;
        for i = setdiff(I, X)
            Y = [X; i];                        
            
            Z = -1;
            if w(i) ~= 0 && chunks(i) > 0
                Z = ((double(chunks(i)) / 1e9 ) / sum(double(chunks(Y)) / 1e9)) *  (1 / (n(i) * w(i)));
            end
            
            if Z > Z_old
                best = i;
                Z_old = Z;
            end
        end
        
        if best ~= 0
            X = [X; best];
        end
        
        K_old = K; 
        B = sum(double(chunks(X)) / 1e9);
        SX = sum(S(X));        
        a = ((double(chunks(X)) ./ 1e9) ./ B);
        b = 1 ./ (w(X) .* n(X));
        c = exp(-(SX-1100)^2 / (2 * 1e3));
        K = sum( a .* b ) * c;                  
    end

    X = setdiff(X, best); %remove last entry since caused a diminishing of K
    
    tmp = zeros(N,1);
    tmp(X) = 1;
    X = tmp;
    
    
    %% for debug, plot BS disposition
    if DEBUG
        figure;
        hold on;
        grid on;
        for i = 1:N
            plot(allBS{i}.pos(1), allBS{i}.pos(2), '*')
            %l = abs(allBS{i}.pos(2) - allBS{i}.sharedData.UE.pos(2)) * sin(pi/3) / sin(pi/6);
            l = S(i)/2;
            plot([allBS{i}.pos(1) - l , allBS{i}.pos(1)+S(i)-l], [allBS{i}.pos(2), allBS{i}.pos(2)] - allBS{i}.ID / 10)
            text(allBS{i}.pos(1), allBS{i}.pos(2)+1, strcat('ID: ', int2str(allBS{i}.ID)));
            text(allBS{i}.pos(1), allBS{i}.pos(2)+2, strcat('w: ', num2str(n(i)*w(i)))); %weigth of BS
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