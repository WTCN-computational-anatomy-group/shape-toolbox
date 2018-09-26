function answer = loopK(N, K)

    time_K = K*( (N+1)/128 + 1/120);
    time_N = N*( (K+1)/128 + 1/120);
    fprintf('loop over K: %f\nloop over N: %f\n', time_K, time_N);
    
    answer = time_K < time_N;

end