function saveits_pp_lay(dir1,file,use_ICCG,use_DICCG,use_POD,dpod,ts,dv,preport,last, varargin)
opt   = struct('ex', [], ...
               'per', []);
[opt] = merge_options(opt, varargin{:});
ex    = opt.ex;
per   = opt.per;

lp = length(dpod);
text = [dir1 file];

for i=1:ts
    iter(i,1)=preport(1,i).iter;
end

if ~exist(text, 'file')
    fileID = fopen(text,'w');
    fprintf(fileID,'\\begin{table}\\small\\centering\n');
    fprintf(fileID,'\\caption{Number of iterations.}\\label{table:it_}   \n');
    fprintf(fileID,' \\begin{tabular}{lllllll}\n');
    fprintf(fileID,'\\hline\\noalign{\\smallskip} \n');
    fprintf(fileID,' $\\frac{\\kappa_2}{\\kappa_1}$ &d &Total& \\multicolumn{2}{c}{DICCG} & Total & \\%% of \\\\ \n');
    fprintf(fileID,'       &   & ICCG     &  ICCG&DICCG &DICCG& ICCG \\\\ \n');
    fprintf(fileID,'\\noalign{\\smallskip}\\hline\\noalign{\\smallskip}\n');
    fclose(fileID);
    
end

if (use_ICCG)
    fileID = fopen(text,'a');
    ttits = sum(iter);
    save([dir1  'ttits.mat'],'ttits')
    fprintf(fileID,['$10^{' num2str(per) '} $&' num2str(dv) '&' num2str(ttits) '&' num2str(ttits) '& 0 ' '&' num2str(ttits)  '&' num2str(1) '\\\\ \n']);
    
else if (use_DICCG)
        if ex
            
            fileID = fopen(text,'a');
            fprintf(fileID,'\\noalign{\\smallskip}\\hline\\noalign{\\smallskip}\n');
            fprintf(fileID,['\\multicolumn{7}{c}{' ex '}\\\\ \n']);
            fprintf(fileID,'\\noalign{\\smallskip}\\hline\\noalign{\\smallskip}\n');
            fclose(fileID);
        end
        icits=sum(iter(1:dv));
        dicits=sum(iter(dv+1:ts));
        tits = sum(iter);
        load([dir1   'ttits.mat'],'ttits')
        perc=round(mean(tits*100/ttits));
        if (use_POD)
            if ex
                fileID = fopen(text,'a');
                fprintf(fileID,'\\noalign{\\smallskip}\\hline\\noalign{\\smallskip}\n');
               fprintf(fileID,['\\multicolumn{7}{c}{' ' POD  deflation vectors}\\\\ \n']);
                fprintf(fileID,'\\noalign{\\smallskip}\\hline\\noalign{\\smallskip}\n');
                fclose(fileID);
            end
            fileID = fopen(text,'a');
            fprintf(fileID,'\\hline  \n');
            fprintf(fileID,['$10^{' num2str(per) '} $&' num2str(numel(dpod)) '&' num2str(ttits)  '&' num2str(icits) '&' num2str(dicits) '&' num2str(tits)  '&' num2str(perc) '\\\\ \n']);
            fprintf(fileID,'\\hline\\noalign{\\smallskip} \n');
            fprintf(fileID,'\\end{tabular} \n');
            fprintf(fileID,'\\end{table}  \n');
            fclose(fileID);
            
            
        else
  %          if ex
  %              fileID = fopen(text,'a');
  %               fprintf(fileID,'\\noalign{\\smallskip}\\hline\\noalign{\\smallskip}\n');
  %              fprintf(fileID,'\\multicolumn{7}{c}{' ex '}\\\\ \n');
  %              fprintf(fileID,'\\noalign{\\smallskip}\\hline\\noalign{\\smallskip}\n');
  %              fclose(fileID);
 %           end
            fileID = fopen(text,'a');
            fprintf(fileID,'\\hline \n');
            fprintf(fileID,['$10^{' num2str(per) '} $&' num2str(dv) '&' num2str(ttits)  '&' num2str(icits) '&' num2str(dicits) '&' num2str(tits)  '&' num2str(perc) '\\\\ \n']);
            fclose(fileID);
            
        end
        if (last)
            fileID = fopen(text,'a');
            fprintf(fileID,'\\hline\\noalign{\\smallskip} \n');
            fprintf(fileID,'\\end{tabular} \n');
            fprintf(fileID,'\\end{table}  \n');
            fclose(fileID);
        end
        
    end
end


end

