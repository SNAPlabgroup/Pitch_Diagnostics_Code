function [x_out, sig_rms] = make_F0DL_stim_ABAB(f0_A, f0_B, rank, nharms, phi, noise_on, dur, ramp, isi, fs, db_main, db_flank, sigrms, rove_A, rove_B)

 
buffer = zeros(round(isi*fs),1);
    
[tone_A, sig_rms_A] = make_F0DL_stim(f0_A, dur, fs, db_main, db_flank, sigrms, rank, nharms, ramp, phi, noise_on, rove_A);
[tone_B, sig_rms_B] = make_F0DL_stim(f0_B, dur, fs, db_main, db_flank, sigrms, rank, nharms, ramp, phi, noise_on, rove_B);

x_out = [tone_A, buffer', tone_B];
x_out = [x_out, buffer', x_out];

sig_rms = [sig_rms_A,sig_rms_B];

end


