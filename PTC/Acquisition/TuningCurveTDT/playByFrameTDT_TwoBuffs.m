%% Just a script to test playing a long array frame-by-frame without
% discontinuities.

% Making up a stimulus of pulsing tones
y = [];
for k = 1:10
    y = [y, maketone(0, 0, 4000, 0.5, 0.02, 48828.125, 0.2)']; %#ok<AGROW>
end
frameDur_samps = 2048;

Lref = 114; % dB SPL output of RZ6-ER2 chain for full-scale MATLAB tone
L = 70;
chanselect = 2; % 1 for left, 2 for right


Nframes = floor(numel(y)/frameDur_samps);

try

    %% Initialize TDT RZ6 for blocked (frame-by-frame) use.
    fig_num=99;
    USB_ch=1;
    FS_tag = 3;
    Fs = 48828.125;
    [f1RZ,RZ,~]=load_play_circuit(FS_tag,fig_num,USB_ch);
    invoke(RZ, 'SetTagVal', 'blocksize',frameDur_samps);
    frameDur = frameDur_samps / Fs;

    % Set desired channel and levels
    if chanselect == 1
        dropA = Lref - L;
        dropB = 120;
        invoke(RZ, 'SetTagVal', 'attA', dropA);
        invoke(RZ, 'SetTagVal', 'attB', dropB);
        invoke(RZ, 'SetTagVal', 'chanselect', 1);
    else
        dropA = 120;
        dropB = Lref - L;
        invoke(RZ, 'SetTagVal', 'attA', dropA);
        invoke(RZ, 'SetTagVal', 'attB', dropB);
        invoke(RZ, 'SetTagVal', 'chanselect', 2);
    end
    %% Play audio out frame by frame.
    t0 = tic;
    count1 = 0;
    count2 = 0;

    for k = 1:Nframes
        idx = (k - 1)*frameDur_samps + (1:frameDur_samps);
        buff = y(idx);
        if mod(k, 2) == 1
            num1 = invoke(RZ, 'GetTagVal', 'numblocks1');
            while num1 < count1
                WaitSecs(0.5*frameDur);
                num1 = invoke(RZ, 'GetTagVal', 'numblocks1');
            end
            count1 = count1 + 1;
            invoke(RZ,'ZeroTag','datain1');
            invoke(RZ, 'WriteTagVEX', 'datain1', 0, 'F32', buff);
        else
            num2 = invoke(RZ, 'GetTagVal', 'numblocks2');
            while num2 < count2
                WaitSecs(0.5*frameDur);
                num2 = invoke(RZ, 'GetTagVal', 'numblocks2');
            end
            count2 = count2 + 1;
            invoke(RZ,'ZeroTag','datain2');
            invoke(RZ, 'WriteTagVEX', 'datain2', 0, 'F32', buff);
        end
        if k == 1
            invoke(RZ, 'SoftTrg', 1); % Start the first frame
        end
    end
    toc(t0);

    %Clearing I/O memory buffers:
    invoke(RZ,'ZeroTag','datain1');
    invoke(RZ,'ZeroTag','datain2');
    close_play_circuit(f1RZ,RZ);

catch me
    %Clearing I/O memory buffers:
    invoke(RZ,'ZeroTag','datain1');
    invoke(RZ,'ZeroTag','datain2');
    close_play_circuit(f1RZ,RZ);
    rethrow(me);
end