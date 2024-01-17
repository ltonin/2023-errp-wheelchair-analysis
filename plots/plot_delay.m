function plot_delay(data,pos, relativ_err)
    %PLOT_DELAY Summary of this function goes here
    %   Detailed explanation goes here
    data = abs(data);
    figure()
    plot(data)
    hold on
    stem(pos, data(pos), 'r')
    stem(pos + relativ_err, data(pos + relativ_err), 'b')
end

