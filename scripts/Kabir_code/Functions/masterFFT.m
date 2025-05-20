function [P1_out,f,complexP1,Phase] = masterFFT(data,options)
    
        %options.Fs = 60; sampling frequency
        %options.padding = 1/0
        %options.padsize = 16/32/64
        %options.hamming = 1/0
    biggness = size(data);
    if length(size(data))==2 && (biggness(1)>1)

        L = size(data,1);             % Length of signal    

        for trl = 1:size(data,2)
            
            dat = data(:,trl)';
            if options.hamming
                dat = (dat'.*hamming(length(dat)))';
            end

            if options.padding
                lpad = (2^nextpow2(L+options.padsize));
                f = options.Fs*(0:(lpad/2))/lpad;
                Y = fft(dat,lpad);
                P2 = abs(Y/lpad).^2;
                P1 = P2(1:lpad/2+1);
            else
                f = options.Fs*(0:(L/2))/L;
                Y = fft(dat);
                P2 = abs(Y/L).^2;
                P1 = P2(1:L/2+1);
            end

            complexP1(:,trl) = Y;
            Phase(:,trl) = angle(Y);

            P1_out(trl,:) = P1;
            P1_out(trl,2:end-1) = 2*P1(2:end-1);
        end
    elseif length(size(data))==3 && (biggness(1)>1)
        for y = 1:size(data,3)
            for yy = 1:size(data,3)
                data_temp = data(:,y,yy);
                L = length(data_temp);             % Length of signal    

                if options.hamming
                    data_temp = (data_temp'.*hamming(length(data_temp)))';
                end

                if options.padding
                    lpad = (2^nextpow2(L+options.padsize));
                    f = options.Fs*(0:(lpad/2))/lpad;
                    Y = fft(data_temp,lpad);
                    P2 = abs(Y/lpad).^2;
                    P1 = P2(1:lpad/2+1);
                else
                    f = options.Fs*(0:(L/2))/L;
                    Y = fft(data_temp);
                    P2 = abs(Y/L).^2;
                    P1 = P2(1:L/2+1);
                end

                complexP1 = Y;
                Phase = angle(Y);

                P1(2:end-1) = 2*P1(2:end-1);
                Pend(y,yy,:) = P1; 
            end
        end
        P1 = Pend;
        P1_out = P1;
    else
        data_temp = data;
        L = length(data_temp);             % Length of signal    

        if options.hamming
            data_temp = (data_temp'.*hamming(length(data_temp)))';
        end

        if options.padding
            lpad = (2^nextpow2(L+options.padsize));
            f = options.Fs*(0:(lpad/2))/lpad;
            Y = fft(data_temp,lpad);
            P2 = abs(Y/lpad).^2.^2;
            P1 = P2(1:lpad/2+1);
        else
            f = options.Fs*(0:(L/2))/L;
            Y = fft(data_temp);
            P2 = abs(Y/L).^2;
            P1 = P2(1:L/2+1);
        end

        complexP1 = Y;
        Phase = angle(Y);

        P1(2:end-1) = 2*P1(2:end-1);
        P1_out = P1;
    end