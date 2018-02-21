%% Create a directory with the variables of the problem to save results

%dir = '/home/wagm/cortes/Localdisk/Research/Report/01_12/';
%dir = '/run/media/taurlin/February/Results/';
dir = '../Results_20/1/';
%dir = '/mnt/sda2/cortes/Research/articles/JCP_18/Results/';
%dir = '/mnt/sda2/cortes/Research/2018/02_2018/Results/';
if(use_wells)
    if(model_SPE)
        folder=[ 'SPE10_' num2str(numel(layers))  'DT_' num2str(DT/day), ...
            'step_' num2str(nstep) 'P_' num2str(P) 'tol_e-' num2str(tol_v) ];
    else
        folder=[ 'per_' num2str(per) 'sz_' num2str(nx)  'layers_', ...
            num2str(nz)  'DT_' num2str(DT/day) 'step_' num2str(nstep), ...
            'P_' num2str(P) 'tol_e-' num2str(tol_v) ];
    end
else
    if(model_SPE)
        folder=[ 'SPE10_' num2str(numel(layers))  'DT_' num2str(DT/day), ...
            'step_' num2str(nstep) 'bc_' num2str(P_b)  'tol_e-' num2str(tol_v)];
    else
        folder=[ 'per_' num2str(per) 'sz_' num2str(nx) '_layers_', ...
            num2str(nz)  'DT_' num2str(DT/day) 'step_' num2str(nstep), ...
            'bc_' num2str(P_b) 'tol_e-' num2str(tol_v)];
    end
end
mkdir([dir], folder)
dir1 = [dir folder '/'];

if use_ICCG
    folder=['ICCG' ];
    if training
        folder=['ICCG_t' ];
    end
else
    if(use_wells)
    folder=['DICCG_dv_' num2str(dv) 'pod_' num2str(numel(dpod)) 'Pp_', ...
        num2str(P_s)];
    else
    folder=['DICCG_dv_' num2str(dv) 'pod_' num2str(numel(dpod)) 'Pb_', ...
        num2str(P_b)];    
    end
end
mkdir([dir1], folder)
dir2 = [dir1 folder '/'];
