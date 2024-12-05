%clear
rr=1
commandwindow
while rr==1
    [keydown,respT] = prm.trigger.btsi.getResponse;
    if keydown>0
        rr=0;
    end

end
fprintf(keydown)