function delays = errp_compute_delay(data, position, compensation)
    delays = ones(length(position),1);

    data = abs(data);
    
    delay_vel = load('mean_delay.mat').delay_vel * compensation;

    for i = 1:length(position)
        j = 0;
        found = false;
        while not(found)
            j = j +1;
            if (data(position(i)+j) > delay_vel)
                found = true;
            end
        end
        delays(i) = j;
    end
    
end