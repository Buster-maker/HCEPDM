function plot_initial_population(pop, V, decision_variables)
    % 决策变量的维度V必须大于或等于2，以便我们绘制前两个决策变量
    if V < 2
        error('V must be greater than or equal to 2 for a 2D plot.');
    end

    % 提取前两个决策变量
    decision_var_x = decision_variables(:, 1);
    decision_var_y = decision_variables(:, 2);

    % 绘制散点图
    figure;
    scatter(decision_var_x, decision_var_y, 'filled');
    title(' Population Distribution');
    xlabel('Decision Variable 1');
    ylabel('Decision Variable 2');
    grid on;
end