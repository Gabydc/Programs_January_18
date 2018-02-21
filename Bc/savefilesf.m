
                if cp==1
                    
                    filecp='cp';
                    B=[dir2   filecp  '.fig'];
                    saveas(hcp,B)
                    B=[dir2   filecp   '.jpg'];
                    saveas(hcp,B)
                end

filetx = ['results.txt'];
%if solver == [1 2 3 4]
%
                saveres(dir1,filetx,def,pod,dpod,per,ts,dv,preport)

%end
                close all
                filews=['workspace'];
                filename=[dir2 filews];
                save(filename)