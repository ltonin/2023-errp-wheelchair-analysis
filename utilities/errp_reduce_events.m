function new_events = errp_reduce_events(events, evtidx, persistance)
    POS = events.POS(evtidx);
    new_POS = POS;

    for i = 1:persistance
        new_POS = setdiff(new_POS, POS + i);
    end

    dif = setdiff(POS, new_POS);
    new_events = evtidx;

    for fIx = 1:length(dif)
        new_events(events.POS == dif(fIx)) = 0;
    end
    new_events = logical(new_events);

end