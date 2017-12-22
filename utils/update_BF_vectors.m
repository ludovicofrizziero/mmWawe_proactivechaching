function [BF_vector_tx, BF_vector_rx] = update_BF_vectors(AoA, AoD, n_tx, n_rx, n_BS_all, serving_BS_idx)
    BF_vector_tx = cell(n_BS_all,1);
    rnd_steering_angle = 2*pi*rand(n_BS_all,1); % random steering angles
    rnd_steering_angle(serving_BS_idx) = AoA;
    
    for bs_idx = 1:n_BS_all
        BF_vector_tx{bs_idx} = compute_BF_vector(n_tx, rnd_steering_angle(bs_idx));
    end

    BF_vector_rx = compute_BF_vector(n_rx,AoD);