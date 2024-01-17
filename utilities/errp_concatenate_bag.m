function [pose, twist, cmdvel, bagevents, baglabels] = errp_concatenate_bag(bagfiles)
    nfiles = length(bagfiles);

    posesize = get_data_size(bagfiles, 'pose', 'x');
    nsample  = sum(posesize(1,:));

    pose.x = nan(nsample,1);
    pose.y = nan(nsample,1);
    pose.z = nan(nsample,1);

    twist.x = nan(nsample,1);
    twist.y = nan(nsample,1);
    twist.z = nan(nsample,1);

    cmdvel.x = nan(nsample,1);
    cmdvel.y = nan(nsample,1);
    cmdvel.z = nan(nsample,1);

    Rk = nan(nsample, 1);                    % Run
   
    TYP = []; POS = []; ERR = [];

    fileseek = 1;
    for fId = 1:nfiles
    
        cfile = bagfiles{fId};
        util_disp_progress(fId, nfiles, '        ');
        
        cdata   = load(cfile);

        % Get current position 
        cstart   = fileseek;
        cstop    = cstart + length(cdata.cmdvel.x) - 1;
        
        % Concatenate the data
        pose.x(cstart:cstop, :) = cdata.pose.x;
        pose.y(cstart:cstop, :) = cdata.pose.y;
        pose.z(cstart:cstop, :) = cdata.pose.z;

        twist.x(cstart:cstop, :) = cdata.twist.x;
        twist.y(cstart:cstop, :) = cdata.twist.y;
        twist.z(cstart:cstop, :) = cdata.twist.z;

        cmdvel.x(cstart:cstop, :) = cdata.cmdvel.x;
        cmdvel.y(cstart:cstop, :) = cdata.cmdvel.y;
        cmdvel.z(cstart:cstop, :) = cdata.cmdvel.z;
    
        TYP = cat(1, TYP, cdata.events.TYP);
        ERR = cat(1, ERR, cdata.events.ERR);
        POS = cat(1, POS, cdata.events.POS + fileseek - 1);

        % Create labels
        Rk(cstart:cstop) = fId;

        % Close the positions
        fileseek = cstop + 1;
    end
    
    bagevents.TYP = TYP;
    bagevents.ERR = ERR;
    bagevents.POS = POS;

    baglabels.samples.Rk = Rk;
    
end


function dsizes = get_data_size(filepaths, variable, subv)

    nfiles = length(filepaths);
    ndimensions = 2;
    dsizes = zeros(ndimensions, nfiles);
    
    for fId = 1:nfiles
        cfilepath = filepaths{fId};
        cdata   = load(cfilepath, variable);
        dsizes(:, fId) = size(cdata.(variable).(subv));
        %eval(['size(cdata.' variable '.' subv ')']);    
    end

end


