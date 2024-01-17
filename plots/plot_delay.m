function plot_delay(data,pos, relativ_err, need_abs)
    % PLOT_DELAY Summary of this function goes here
    %   Detailed explanation goes here
    if need_abs
        data = abs(data);
    end

    figure()
    plot(data)
    hold on
    stem(pos, data(pos), 'r')
    stem(pos + relativ_err, data(pos + relativ_err), 'b')
    
end

