function savefigures(f, title, dir)

B=[dir title '.fig'];
savefig(f,B)
B=[dir title   '.jpg'];
saveas(f,B)
end