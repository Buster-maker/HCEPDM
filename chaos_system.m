% Henon chaotic system
function seeds = chaos_system(pop, V)
    % Henon map parameters
    a = 1.4;
    b = 0.3;
    seeds = zeros(pop, V);
    for i = 1:pop
        x = rand(1, V);
        y = rand(1, V);
        for j = 1:V
            x_new = 1 - a * x(j)^2 + y(j);
            y_new = b * x(j);
            x(j) = x_new;
            y(j) = y_new;
        end
        seeds(i, :) = x;
    end
end