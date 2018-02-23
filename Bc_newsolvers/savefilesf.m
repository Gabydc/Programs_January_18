 use_Cp =0;
                if use_Cp
                    
                    filecp='cp';
                    B=[dir2   filecp  '.fig'];
                    saveas(hcp,B)
                    B=[dir2   filecp   '.jpg'];
                    saveas(hcp,B)
                end

filetx = ['results.txt'];
%if solver == [1 2 3 4]
use_DICCG =def;
use_POD = pod;
                saveres(dir1,filetx,use_DICCG,use_POD,dpod,per,ts,dv,preport)

%end
                close all
                filews=['workspace'];
                filename=[dir2 filews];
                save(filename)