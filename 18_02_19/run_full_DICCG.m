use_DICCG = false;    

for dv = [10 30]
    clearvars -except use_DICCG dv
for use_ICCG = [true, false]
    if (~use_ICCG)
        use_DICCG = true;    
    end
        win_2ph_bc_lay       
    
end
end