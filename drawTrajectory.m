function drawTrajectory(X,t)
%DRAWTRAJECTORY plot each row of X
    % Parse arguments
    [m,n]=size(X);
    if nargin==1
        t=1:n;
    end
    cm = colormap('Lines');
    % hold on
    cla  % clear the figure
    for i=1:m
        % figure;
        hold on;
        plot(t,X(i,:),  'Color', cm(i,:), ...
            'LineWidth', 2);
        plot(t(end),X(i,end), '*k', 'MarkerSize', 10);
    end
    hold off;
end

