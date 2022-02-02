function A = generate_flow_network(n_vertices,p)
% This function generates a random flow network, i.e. a weighted directed  
% graph, that is to be interpreted as a chemical reaction. The parameter p
% is the probability with which an edge is added.

% Initialize adjacency matrix.
A = zeros(n_vertices,n_vertices);

% W.l.o.g. assume that there is a chain of transformations in the order
% 1:n_vertices, with v_1 being the educt and v_n being the product of the
% reaction. Add edges to all vertices besided v_n and specify weights s.t. 
% each vertex has more inflow than outflow.
weights = sort(unifrnd(0.25,1,[n_vertices-2,1]),'descend');
for i = 2:n_vertices-1
    A(i-1,i) = weights(i-1);
end

% Add edges to product with probability p.
edges = rand(n_vertices-1,1) <= p;
if nnz(edges) == 0
    % Add one edge from a random vertex in the case that no edge has been
    % randomly generated above.
    edges(randi(n_vertices-1)) = 1;
end

% In the case that v_{n-1} has not been connected to v_n above, we have to
% add an edge from v_{n-1} to another random vertex. Also make sure that
% inflow > outflow.
if edges(n_vertices-1) == 0
    inflow = sum(A(:,n_vertices-1));
    v = randi(n_vertices-2);
    A(n_vertices-1,v) = unifrnd(0,inflow);
end
    
    
    
% Add weights to the generated edges, making sure that inflow > outflow for
% all vertices besides v_1 and outflow > inflow for v_1.
weights = zeros(nnz(edges),1);
vertices = find(edges);
for i = 1:nnz(edges)
    v = vertices(i);
    if v ~= 1
        inflow = sum(A(:,v));
        outflow = sum(A(v,:));
        dif = inflow - outflow;
        weights(i) = unifrnd(0,dif);
    else
        weights(i) = unifrnd(0,0.25);
    end
end
    
A(edges,n_vertices) = weights;

% Now that all species are connected correctly, the rest of random edges is
% added with probability p. Already existing edges and self-loops are 
% excluded.
existing_edges = A > 0;
rand_edges = rand(n_vertices,n_vertices) <= p;
new_edges = ~existing_edges & rand_edges;
new_edges = tril(new_edges,-1) + triu(new_edges,1)';

% Delete outgoing edges from product.
new_edges(n_vertices,:) = zeros(1,n_vertices);

% Set weights for v_1. Begin with outgoing edges.
out_edges = find(new_edges(1,:));
weights = unifrnd(0,0.25,[1,length(out_edges)]);
A(1,out_edges) = weights;

% Set weights for incoming edges of v_1 with outflow > inflow.
outflow = sum(A(1,:));
in_edges = find(new_edges(:,1));
weights = zeros(length(in_edges),1);
for i = 1:length(in_edges)
    v = in_edges(i);
    % Determine inflow - outflow of connected vertex. The weight of the
    % edge has to be smaller than the outflow of v_1 and the gap between
    % inflow and outflow of the connected node.
    dif = sum(A(:,v)) - sum(A(v,:));
    weights(i) = unifrnd(0,min(outflow,dif));
    outflow = outflow - weights(i);
end
A(in_edges,1) = weights;

% Set edges to zero as they have already been processed.
new_edges(1,:) = zeros(1,n_vertices);
new_edges(:,1) = zeros(n_vertices,1);

% Set weights for intermediate nodes. Make sure that inflow > outflow.
for k = 2:n_vertices
    
    % Set weigths for incoming edges.
    in_edges = find(new_edges(:,k));
    weights = zeros(length(in_edges),1);
    for i = 1:length(in_edges)
        v = in_edges(i);
        % Determine inflow - outflow of connected vertex. The weight of the
        % edge has to be smaller than the gap between inflow and outflow of
        % the connected node.
        dif = sum(A(:,v)) - sum(A(v,:));
        weights(i) = unifrnd(0,dif);
    end
    A(in_edges,k) = weights;
    
    % Set edges that have already been processed to zero.
    new_edges(:,k) = zeros(n_vertices,1);
    
    % Set weigths for outgoing edges. The gap between between inflow and
    % outflow of the current vertex is the upper bound for the total amount
    % of weight that still can be distributed.
    dif = sum(A(:,k)) - sum(A(k,:));
    out_edges = find(new_edges(k,:));
    weights = zeros(1,length(out_edges));
    for i = 1:length(out_edges)
        weights(i) = unifrnd(0,dif);
        dif = dif - weights(i);
    end
    A(k,out_edges) = weights;
    new_edges(k,:) = zeros(1,n_vertices);
end

end