
%% Define Solver
training   = false;


% use_ICCG   = false;
% use_DICCG  = true;
% use_POD    = false;
if (use_POD)
    last = true;  % Last experiment, to close the table
else
    last = false;
end
plot_sol   = false;
save_res   = true;
use_wells  = false;
model_SPE  = false;
use_cp     = false;
window     = true;
use_agmg   = false;


if use_DICCG
    training = false;
end
if use_ICCG
    use_POD = false;
end