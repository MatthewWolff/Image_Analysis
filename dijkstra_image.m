function [dist] = dijkstra_image(map, source, goal)
%assumes that source and goal are (column, row) coordinates, not linear indices
% Based off this Java code: http://cs.fit.edu/~ryan/java/programs/graph/Dijkstra-java.html

assert(length(size(map)) == 2)
get_weight = @(c) map(c(1),c(2));
find_neighbors = @(x,y) [x-1,y-1;x-1,y;x-1,y+1;x,y-1;x,y+1;x+1,y+1;x+1,y;x+1,y-1];
nodes = find(map(:) ~= 0); % linear index of all nodes
node = @(location) ismember(nodes, location,'rows'); % uses node linear index to find position in nodes
coord2lin = @(x) sub2ind(size(map), x(1), x(2));

origin = coord2lin(flip(source));
terminal = coord2lin(flip(goal));

visited = zeros(length(nodes),1); % all false
distance = Inf(length(nodes),1);
predecessor = zeros(length(nodes),1);
distance(node(origin)) = 0;




for i = 1:length(nodes)
    tempA = Inf;
    tempB = -1;   % graph not connected or no unvisited vertices
    for j = 1:length(nodes)
        if (~visited(j) && distance(j) < tempA) %if unvisited
            tempB=nodes(j); % remember closest node?
            tempA=distance(j); % remembers closest node distance?
        end
    end
    next = tempB;
    visited(node(next)) = true;
    
    [x,y] = ind2sub(size(map,1), next); %convert index to coords for neighbor finding
    neighbors = find_neighbors(x,y);
    neighbors = neighbors(neighbors(:,1) > 0 & neighbors(:,2) > 0, :); %remove impossible neighbors
    for j = 1 : length(neighbors)    % look thru neighbors
        n = neighbors(j,:); % this iteration's neighbor
        local_dist = distance(node(next)) + get_weight(n);%distance from point to neighbor == neighbor value
        if (distance(node(coord2lin(n))) > local_dist)
            distance(node(coord2lin(n))) = local_dist;
            predecessor(node(coord2lin(n))) = next;
        end
    end
end

% DISPLAY PATH
% path = [];
% trace = terminal;
% while (trace ~= origin) 
%     path =  [path, predecessor(node(trace))]; % assembles backwards
%     trace = predecessor(node(trace));
% end
% display(flip(path)
dist = distance(node(terminal));

end