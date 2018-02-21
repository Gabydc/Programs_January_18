function saveits(dir1,file,use_ICCG,use_DICCG,use_POD,dpod,ts,dv,preport,last)

text = [dir1 file];
for i=1:ts
    its(i,1)=preport(1,i).iter;
end

if ~exist(text, 'file')
    fileID = fopen(text,'w');
    fprintf(fileID,'\\begin{table}[!ht]\\centering\n');
    fprintf(fileID,'\\begin{minipage}{0.5\\textwidth}\n');
    fprintf(fileID,' \\centering\n');
    fprintf(fileID,'\\begin{tabular}{ ||c|c|c|c|c||} \n');
    fprintf(fileID,'\\hline\n');
    fprintf(fileID,'Total&Method & ICCG&DICCG &\\%% of ICCG\\\\ \n');
    fprintf(fileID,'ICCG& & &Iterations\\\\ \n');
    fclose(fileID);
    
end

if (use_ICCG)
    fileID = fopen(text,'a');
    ttits = sum(its);
    save([dir1  'ttits.mat'],'ttits')
   fprintf(fileID,[ num2str(ttits) '& ICCG&' num2str(ttits)  '&' num2str(1) '\\\\ \n']);

else if (use_DICCG)
        
        tits = sum(its);
        load([dir1   'ttits.mat'],'ttits')
        perc=round(mean(tits*100/ttits));
        if (use_POD)
            fileID = fopen(text,'a');
            fprintf(fileID,'\\hline  \n');
            fprintf(fileID,[ num2str(ttits) '& DICCG$_{POD_{' num2str(numel(dpod)) '}}$&' num2str(tits)  '&' num2str(perc) ' \\\\ \n']);
            fclose(fileID);
            if (last)
                fileID = fopen(text,'a');
                fprintf(fileID,'\\hline  \n');
                fprintf(fileID,'\\end{tabular} \n');
                fprintf(fileID,'\\caption{Comparison between the ICCC and DICCG methods of the average number of linear iterations. }\\label{table:litertot} \n');
                fprintf(fileID,'\\end{minipage}  \n');
                fprintf(fileID,'\\end{table}  \n');
                fclose(fileID);
            end
            
        else
            fileID = fopen(text,'a');
            fprintf(fileID,'\\hline \n');
            fprintf(fileID,[ num2str(ttits) '& DICCG$_{' num2str(dv) '}$&' num2str(tits)  '&' num2str(perc) '\\\\ \n']);
            fclose(fileID);
            
        end
        
    end
end


end

