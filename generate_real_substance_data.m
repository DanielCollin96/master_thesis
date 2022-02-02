function [A, W, H, K, tspan, freq] = generate_real_substance_data()
% This function generates a Raman spectral data set example, using real
% substances as species.

% Reaction coefficients.
K = [...
    -0.53,  0.02,  0.0,   0,  0.0;
     0.53, -0.66,  0.25,  0,  0.0;
     0.0,   0.43, -0.36,  0,  0.1;
     0.0,   0.21,  0.0,   0,  0.0;
     0.0,   0.0,   0.11,  0, -0.1];

% Time interval.
tspan = linspace(0,50,151);

% Initial concentrations.
y0 = [1,0,0,0,0]';

% Compute relative concentrations to obtain H.
odefun = @(t,y) K*y;
[~,H] = ode45(odefun, tspan, y0);
H = H';

% Construct W from the original data.
spec = load('data\trans_dichlorethen.dat');
W(:,1) = spec(:,2);
spec = load('data\toluene.dat');
W(:,2) = spec(:,2);
spec = load('data\diethylether.dat');
W(:,3) = spec(:,2);
spec = load('data\cis_dichlorethen.dat');
W(:,4) = spec(:,2);
spec = load('data\acetonitrile.dat');
W(:,5) = spec(:,2);

% Store frequencies.
freq = spec(:,1);

clear spec;

% Compute A.
A = W*H;

end