function ellipse = cov_ellipse(mu, P, scale, n)
    % Parameter used to define ellispe
    t = linspace(0, 2*pi, n);

    % Base circle to be transformed
    ellipse = [cos(t); sin(t)];
    
    % Getting Eigenvalues and Eigenvectors of covariance
    [V, D] = eig(P);

    ellipse = V * scale * sqrt(D) * ellipse + mu;

end