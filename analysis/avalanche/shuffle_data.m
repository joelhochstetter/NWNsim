function data_randomly_placed =  shuffle_data(data)

     data_randomly_placed = zeros(size(data));

    rand_pos = randperm(numel(data)); %array of random positions
    % new array with original data randomly distributed 
    for k = 1:numel(data)
        data_randomly_placed(k) = data(rand_pos(k));
    end
    
end