function [K,h0] = generate_kinetics(n_species,p)
% This function generates a random reaction kinetics by creating a flow
% network and transforming it into a reaction coefficient matrix. The
% parameter p is the probability with which an edge is added.

% Get adjacency matrix of a random flow network.
A = generate_flow_network(n_species,p);

% Transpose and calculate diagonal entries.
K = A'; 
for i = 1:n_species
    K(i,i) = -sum(K(:,i));
end

% Scale matrix such that all entries are smaller or equal 1.
if max(max(abs(K))) > 1
    K = K / max(max(abs(K)));
end

% Set species a as only educt of reaction.
h0 = zeros(n_species,1);
h0(1) = 1;

end
