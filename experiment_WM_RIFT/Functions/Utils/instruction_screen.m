function instruction_screen(prm,txt)

Screen('glPushMatrix', prm.monitor.window);
Screen('glLoadIdentity',prm.monitor.window);

default('txt',prm.exp.instructions)
Screen('Flip', prm.monitor.window);  
DrawFormattedText(prm.monitor.window, txt, ...
                   'center','center', prm.monitor.white, [], [], [], 1.5);

Screen('Flip', prm.monitor.window);    
Screen('glPopMatrix', prm.monitor.window);

pause(.5);KbWait; 
Screen('Flip', prm.monitor.window);

