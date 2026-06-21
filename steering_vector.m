function a = steering_vector(N, thetha)
    thetha_rad = deg2rad(thetha);
    a = zeros(N,1);
    for m = 1:N
        idx = 2*m - N - 1;
        a(m) = exp(-1j*idx*pi*sin(thetha_rad)/2);
    end
end