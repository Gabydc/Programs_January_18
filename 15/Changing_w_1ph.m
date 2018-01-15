%% Changing wells
%Create a well structure to support multiple report steps
niw = numel(W);
[newW{1:niw+1,1}] = deal(W);
W0 = newW;  % Wells with the same conditions as the original set

% We create a vector with the bhp or rate values for all the timesteps
Value = zeros(numel(W),niw);
%Valr = zeros(numel(W),steps);


for w = 1:numel(newW),
    for  i = 1:4
        if w == i
        W1{w}(1)=newW{w}(w);
        end
       
    end
    
W1{w}(2)=newW{1}(5);    
end


 % New wells properties

clear newW


        