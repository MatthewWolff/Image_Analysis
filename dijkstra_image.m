function [dist, path] = dijkstra_image(map, source, goal)
% assumes that source and goal are (column, row) coordinates, not linear indices
% assumes that there is only 1 source, but 1 or more goals.
% Based off this Java code: http://cs.fit.edu/~ryan/java/programs/graph/Dijkstra-java.html

assert(length(size(map)) == 2)
assert(size(goal,1) > 0)
assert(length(source) == 2 && size(source,1) == 1)

get_weight = @(c) map(c(1),c(2));
find_neighbors = @(x,y) [x-1,y-1;x-1,y;x-1,y+1;x,y-1;x,y+1;x+1,y+1;x+1,y;x+1,y-1];
nodes = find(map); % linear indices of all nodes
node = @(location) ismember(nodes, location,'rows'); % uses node linear index to find position in nodes
coord2lin = @(x) sub2ind(size(map), x(1), x(2));

% change given coordinates into 'proper' format (linear index)
origin = coord2lin(fliplr(source));
if(size(goal,1) > 1)
    terminal = fliplr(goal);
    terminal = sub2ind(size(map), terminal(:,1), terminal(:,2));
else
    terminal = coord2lin(flip(goal));
end

visited = zeros(length(nodes),1); % all false
distance = Inf(length(nodes),1); % all infinite until otherwise known
predecessor = zeros(length(nodes),1); % the list to help us find paths
distance(node(origin)) = 0; % the distance to each node of the node list
        
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
trace = terminal;
if(size(trace,1) == 1) % one goal
    path = [];
else
    path = cell(size(trace)); %multiple goals
end

if(size(trace,1) == 1) % one goal
    while (trace ~= origin)
        path =  [predecessor(node(trace)); path]; % assembles backwards
        trace = predecessor(node(trace));
    end
else % multiple goals
    for i = 1:size(trace,1)
        temp_path = [];
        while (trace(i,:) ~= origin)
            temp_path = [predecessor(node(trace(i,:)));temp_path]; % assembles backwards
            trace(i,:) = predecessor(node(trace(i,:)));
        end
        path{i} = temp_path;
    end
end
dist = distance(ismember(nodes, terminal,'rows'));

end

%% Dijkstra's Method using a Priority Queue (faster, but unsupported)
% MATLAB allows you to call Java code but not all the parts that I need...
% NEED: to supply a comparator so that the priority queue can order itself
% pq = java.util.PriorityQueue;
% pq.add({0,origin});
% distance(node(origin)) = 0; % the distance to each node of the node list
% 
% 
% while(~pq.isEmpty())
%     head = pq.poll();
%     current = head(2);
%     
%     [x,y] = ind2sub(size(map,1), current); %convert index to coords for neighbor finding
%     neighbors = find_neighbors(x,y);
%     neighbors = neighbors(neighbors(:,1) > 0 & neighbors(:,2) > 0, :); %remove impossible neighbors
%     
%     for j = 1 : length(neighbors)    % look thru neighbors
%         n = neighbors(j,:); % this iteration's neighbor
%         local_dist = distance(node(current)) + get_weight(n);%distance from point to neighbor == neighbor value
%         display(distance(node(coord2lin(n))) > local_dist)
%         if (distance(node(coord2lin(n))) > local_dist)
%             distance(node(coord2lin(n))) = local_dist;
%             pq.add({distance(node(coord2lin(n))),coord2lin(n)});
%         end
%     end
% end
% fprintf('Vertex   Distance from Source\n');
% for i = 1:length(nodes)
%     fprintf('%d \t\t %d\n', i, distance(i));
% end