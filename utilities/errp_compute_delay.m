function delays = errp_compute_delay(data, position)
    % Function that return when the data in position vaires more that a
    % treshold
    delays = ones(length(position),1);

    data = abs(data);
    
    f = fft(data(~isnan(data)));
    Fs = 512;
     T =1 /Fs;
     L = length(f);
     t = (0:L-1)*T;
    %figure();
    %plot(Fs/L*(0:L-1),(abs(f)),"LineWidth",3);

    Ft = 0.5;
    Ftn = Ft/(Fs/2);
    n = 500;
    b = fir1(n, Ftn, 'low');
    %data = filter(b,1,data);
   
    %figure()
    %plot(x1)
    %hold on
    %plot(data)

    
    % Figure out when there is the first peak after the command
    for i = 1:length(position)
        found = false;
        j = 0;
        tmp_data = data(position(i));
        
        while not(found)
            j = j+1; 

           
            c_max = max(tmp_data,data(position(i)+j));

            if (tmp_data < c_max)
                tmp_data = c_max;
            else
                second = 0;
                l = 0;
                while (second < 700)
                    second = second + 1;

                    tmp_max = max(tmp_data, data(position(i)+j+second));
                    if ( tmp_max > tmp_data)
                        tmp_data = tmp_max;
                        l = second;
                    end

                end

                found = true;
                delays(i) = j+l;
            end
        end

    end
end

function delays = mtd1()
    delays = ones(length(position),1);
    
    delay_mean = load('mean_delay.mat').delay_mean;

    frame = 1/512;
    pulse = delay_mean/frame;
    delay_mean_pulsation = round(pulse);

    delays = delays*delay_mean_pulsation; 
end

function delays = mtd2()
    delays = ones(length(position),1);

   
    dv = 0.02;

    data = abs(data);

    for i = 1:length(position)
        found = false;
        j = 0;
        while not(found)
            j = j+1; 
            m_min = min(data(position(i)), data(position(i)+j));
            m_max = max(data(position(i)), data(position(i)+j));
            current_error = (m_max-m_min) ^ 2;
            if current_error > dv
                found = true;
                delays(i) = j;
            end
        end

    end
end



% 
%     % Plot where the error is expected dv = 0.1; % how much the velocity
%     is different in order to registrate the alignment
% 
%     ErrorsLx = evtbag.POS(find(evtbag.TYP == prs.ErrorLx)); ErrorsRx =
%     evtbag.POS(find(evtbag.TYP == prs.ErrorRx));
% 
%     NoReleaseLx = evtbag.POS(find(evtbag.TYP == prs.NoReleaseLx));
%     NoReleaseRx = evtbag.POS(find(evtbag.TYP == prs.NoReleaseRx));
% 
%     errorsLx = find_delays(bag, ErrorsLx, dv); errorsRx =
%     find_delays(bag, ErrorsRx, dv);
% 
%     noReleaseLx = find_delays(bag, NoReleaseLx, dv); noReleaseRx =
%     find_delays(bag, NoReleaseRx, dv);
% 
%     errors = [errorsLx; errorsRx];%, noReleaseLx, noReleaseRx]; Errors =
%     [ErrorsLx; ErrorsRx];%; NoReleaseLx; NoReleaseRx];
% 
%     util_bdisp('[vel] + Find alignment between gdf and bag files:');
%     disp(['     | Alignment (LX): ' num2str(mean(errorsLx))]); disp(['
%     | std (LX): ' num2str(std(errorsLx))]); disp(['     | Alignment (RX):
%     ' num2str(mean(errorsRx))]); disp(['     | std (RX): '
%     num2str(std(errorsRx))]); disp(['     | Alignment (Full): '
%     num2str(mean(errors))]); disp(['     | std (full): '
%     num2str(std(errors))]);
% 
%     vel_alignment = round(mean(errors));
% 
%     figure(); plot(bag.twist.z); for i = 1:length(Errors)
%         xline(Errors(i), 'r'); xline(Errors(i)+alignment, 'b');
%         xline(Errors(i)+vel_alignment, 'black');
%         xline(Errors(i)+errors(i), 'g');
%     end
% end
% 
% 
% function alignment = errp_find_alignment_vel_angular(bag, evtbag, evtgdf,
% prs)
% 
%     % Get all events from gdf file events_of_interest =
%     unique(evtgdf.TYP);
% 
%     % Get event index from bag evtbag_index =
%     get_event_index(events_of_interest, evtbag.TYP);
% 
%     % Extract events from bag structure evtbag.TYP =
%     evtbag.TYP(evtbag_index); evtbag.POS = evtbag.POS(evtbag_index);
%     evtbag.ERR = evtbag.ERR(evtbag_index);
% 
%     % Find GDF sequence in bag event list start_idx =
%     find_sequence(evtgdf.TYP, evtbag.TYP);
% 
%     % If sequence has been found, than use the first event to compute the
%     % displacement if(isempty(start_idx) == false)
%         alignment = evtbag.POS(start_idx) - evtgdf.POS(1);
%     end
% 
%     vel_alignment = compute_vel_alignment(bag, evtbag, prs, alignment);
%     alignment = vel_alignment;
% 
% end
% 
% function delays = find_delays(bag, ErrorsPos, dv)
%     delays = zeros(length(ErrorsPos),1);
% 
%     for i = 1:length(ErrorsPos)
%         found = false; j = 0; while not(found)
%             j = j+1; current_error = (bag.twist.z(ErrorsPos(i)) -
%             bag.twist.z(ErrorsPos(i)+j))^2; if current_error > dv
%                 found = true; delays(i) = j;
%             end
%         end
% 
%     end
% end