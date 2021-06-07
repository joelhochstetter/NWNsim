function events = loadmat(eventFolder)
    e = load(eventFolder, 'events');
    events = e.events;
end