function soundTDTR(y, RZ, drop)
% Play sound in the right ear via TDT with optional analog attenuation

stimlength = size(y, 1);

stimTrigger = 1; % Does nothing here
invoke(RZ, 'SetTagVal', 'trigval', stimTrigger);
invoke(RZ, 'SetTagVal', 'nsamps', stimlength);
invoke(RZ, 'WriteTagVEX', 'datainL', 0, 'F32', zeros(size(y))); %write to buffer left ear
invoke(RZ, 'WriteTagVEX', 'datainR', 0, 'F32', y); %write to buffer right ear

invoke(RZ, 'SetTagVal', 'attA', 120); %setting analog attenuation L
invoke(RZ, 'SetTagVal', 'attB', drop); %setting analog attenuation R

WaitSecs(0.05); % Just giving time for data to be written into buffer
%Start playing from the buffer:
invoke(RZ, 'SoftTrg', 1); %Playback trigger