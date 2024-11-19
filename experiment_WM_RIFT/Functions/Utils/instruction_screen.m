function instruction_screen(prm,txt)

Screen('glPushMatrix', prm.w.Number);
Screen('glLoadIdentity',prm.w.Number);

%default('txt',prm.exp.instructions)
Screen('Flip', prm.w.Number);  
DrawFormattedText(prm.w.Number, txt, ...
                   'center','center', prm.monitor.white, [], [], [], 1.5);

Screen('Flip', prm.w.Number);    
Screen('glPopMatrix', prm.w.Number);

warning('In instruction. Press Spacebar to continue.');
waitforspace; waitfornokey;