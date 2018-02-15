%% Changing wells
%Create a well structure to support multiple report steps
[newW{1:nstep,1}] = deal(W);
W0 = newW;  % Wells with the same conditions as the original set

% We create a vector with the bhp or rate values for all the timesteps
Value = zeros(numel(W),nstep);
%Valr = zeros(numel(W),steps);

%Then we'll assign the time-dependent controls
% we select the number of timesteps that have the same pressure value (tch)

csteps = round(nstep/tch);  % To assure that the number of timesteps is integer

% We asign the pressure value to the vectors
h = 0; %counter
while tch*(h+1) < nstep
    h=h+1;
    Value(1,(h-1)*tch+1:tch*h) = rand;
    Value(3,(h-1)*tch+1:tch*h) = rand;
end
Value(2,:) = 1-Value(3,:);
Value(4,:) = 1-Value(1,:);
Value(5,:) = I;
% We can also construct a ramp (increasing values)
%Value(i,1:steps) = linspace(0.5, 0.8, steps);


for w = 1:numel(newW),
    for i = 1:2
        newW{w}(i).type = 'bhp';
        newW{w}(i).val= (P/2*Value(i,w)+P/2)*units_P;
    end
    for i = 3:4
        newW{w}(i).type = 'bhp';
        newW{w}(i).val= (P/2*Value(i,w)+P/2)*units_P;
    end
end


W1=newW; % New wells properties
clear newW


        