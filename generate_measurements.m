function noisy_measurements = generate_measurements(h_cords, R, n)
    noise = mvnrnd(zeros(1,3), R, n);
    noisy_measurements = h_cords + noise;
end