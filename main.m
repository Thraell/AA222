nbIndiv = 100; nbSelec = 25; nbRules = [1;2;3;4;5]; bp = 0.05; maxGen = 20;

for i = 1:1:length(nbRules)

    fprintf('Start :')
    [b(i),Fcb(i)] = gaWBCD(nbIndiv, nbSelec, nbRules(i), bp, maxGen);
    fprintf('End \n')

end