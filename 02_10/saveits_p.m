function saveits_p(dir1,file,use_ICCG,use_DICCG,use_POD,dpod,ts,dv,preport,last, varargin)
opt   = struct('ex', []);
[opt] = merge_options(opt, varargin{:});
ex    = opt.ex;

lp = length(dpod);
text = [dir1 file];

for i=1:ts
    its(i,1)=preport(1,i).iter;
end

if ~exist(text, 'file')
    fileID = fopen(text,'w');
    fprintf(fileID,'\\begin{table}\\small\\centering\n');
    fprintf(fileID,'\\caption{Number of iterations.}\\label{table:it_}   \n');
    fprintf(fileID,' \\begin{tabular}{lllllll}\n');
    fprintf(fileID,'\\hline\\noalign{\\smallskip} \n');
    fprintf(fileID,' $\\frac{\\kappa_2}{\\kappa_1}$ &d &Total& \\multicolumn{2}{c}{DICCG} &Total&\\% of \\\\ \n');
    fprintf(fileID,'       &   & ICCG     &  ICCG&DICCG &DICCG& ICCG\\\\ \n');
    fprintf(fileID,'\\noalign{\\smallskip}\\hline\\noalign{\\smallskip}\n');
    fprintf(fileID,'\\multicolumn{7}{c}{}\\\\ \n');
    fprintf(fileID,'\\noalign{\\smallskip}\\hline\\noalign{\\smallskip}\n');   
    fclose(fileID);
    
end

if (use_ICCG)
    fileID = fopen(text,'a');
    ttits = sum(its);
    save([dir1  'ttits.mat'],'ttits')
   fprintf(fileID,[ num2str(ttits) '& ICCG&' '& &' num2str(ttits)  '&' num2str(1) '\\\\ \n']);

else if (use_DICCG)
        icits=sum(its(1:dv));
        dicits=sum(its(dv+1:ts));
        tits = sum(its);
        load([dir1   'ttits.mat'],'ttits')
        perc=round(mean(tits*100/ttits));
        if (use_POD)
            fileID = fopen(text,'a');
            fprintf(fileID,'\\hline  \n');
            fprintf(fileID,[ num2str(ttits) '& DICCG$_{POD_{' num2str(numel(dpod)) '}}$&' num2str(icits) '&' num2str(dicits) '&' num2str(tits)   '&' num2str(perc) ' \\\\ \n']);
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
            fprintf(fileID,[ num2str(ttits) '& DICCG$_{' num2str(dv) '}$&' num2str(icits) '&' num2str(dicits) '&' num2str(tits)    '&' num2str(perc) '\\\\ \n']);
            fclose(fileID);
            
        end
        
    end
end


end

