function newcell = batch_regexprep(oldcell,oldvals,newvals)
%does a bunch of find-and-replace operations on a cell array using
%regexrep. Give it a cell array of old values {'a','b','c','d'} and a cell
%array of new values to change them to {'x','y','z','q'}. Replacements are
%applied sequentially, so be careful!

newcell = oldcell;

for i = 1:length(oldvals)
    
    newcell = regexprep(newcell,oldvals{i},newvals{i});
    
    
end

end