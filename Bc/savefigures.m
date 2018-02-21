%%If we want to save the files/graphs, it's neccesary to run this part
function savefigures(f, title, dir)    
               
                B=[dir title '.fig'];
                savefig(f,B)
                 B=[dir title   '.jpg'];
                 saveas(f,B)
end