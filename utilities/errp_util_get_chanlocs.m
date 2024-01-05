function chanlocs = errp_util_get_chanlocs(labels, allchanlocs)

    [~, idx] = ismember(lower(labels), lower({allchanlocs.labels}));

    chanlocs = allchanlocs(idx);
    

end