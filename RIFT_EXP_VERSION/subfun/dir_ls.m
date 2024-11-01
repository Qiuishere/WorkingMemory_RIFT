function d = dir_ls(inp)
if nargin==0
    a = dir;
else
    a = dir(inp);
end
for dd = 1:numel(a)
    ok(dd) = 1;
    if strcmp(a(dd).name,'.') || strcmp(a(dd).name,'..')
        ok(dd) = 0;
    end
end
d = a(ok==1);
        
    
