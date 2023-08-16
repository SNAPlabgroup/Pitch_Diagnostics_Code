clear;
clc;

Nfiles = 8; % Number of files per parameter level
minparam = 1;
maxparam = 40;
params = minparam:maxparam;
m = logspace(log10(.25), log10(20), numel(params))/100; %F0DL to test
dur = 0.2;
pad = 0.3;
isi = .5;
fs = 48828.125;
sigrms = 0.020;
db_main = 70;
db_flank = 64;
ramp = 0.02;
% ranks = [4,5,6,7,8,9,10,11,13,16];
ranks = 2:2:12;
nharms = 4;
f0 = 103;
noise_on = 1;

%Alt phase
phi = ones(nharms+2, 1) * pi/2;
phi(1:2:end) = 0;
phi = phi + pi/4; % offset needed to ensure stim starts at zero during alt phase.

ind = 0;
for r = 1:length(ranks)
    rank = ranks(r);
    for k = 1:numel(params)
        difLimen = m(k);
        for nf = 1:Nfiles
            answer = randi(3);
            param = params(k);

            rove = rand()>.5;
            rove_d = rand(2,1)>.5;
            
            ind = ind + 1;
            rove_check(ind,1:3) = [rove,rove_d'];
            choice_check(ind) = answer;
                
            %randomize (no replacement up/down)
            shift_ud = rand()>.5;

            disp(m)
            if shift_ud
                [sig, ~] = make_F0DL_stim_ABAB(f0, f0*(1+m(k)), rank, nharms, phi, noise_on, dur, ramp, pad, fs, db_main, db_flank, sigrms, rove, ~rove);
                for d = 1:2
                    [dummy(d,:), ~] = make_F0DL_stim_ABAB(f0, f0, rank, nharms, phi, noise_on, dur, ramp, pad, fs, db_main, db_flank, sigrms, rove_d(d), ~rove_d(d));
                end            
            else
                [sig, ~] = make_F0DL_stim_ABAB(f0, f0*(1-m(k)), rank, nharms, phi, noise_on, dur, ramp, pad, fs, db_main, db_flank, sigrms, rove_d, ~rove_d);
                for d = 1:2
                    [dummy(d,:), ~] = make_F0DL_stim_ABAB(f0, f0, rank, nharms, phi, noise_on, dur, ramp, pad, fs, db_main, db_flank, sigrms, rove_d(d), ~rove_d(d));
                end
            end
            
            buff = zeros(2,round(isi*fs));
            sig = [sig, buff(1,:)];
            dummy = [dummy,buff];

            switch answer
                case 1
                    y = [sig, dummy(1,:), dummy(2,:)];
                case 2
                    y = [dummy(1,:), sig, dummy(2,:)];
                case 3
                    y = [dummy(1,:), dummy(2,:), sig];
            end

            fname = strcat('./trialaudio/trial',num2str(param),...
                '_',num2str(rank),'_', num2str(nf), '.mat');
            save(fname, 'y', 'fs', 'difLimen', 'rank','phi','nharms','db_main','db_flank', 'param', 'answer');
            
            clear dummy sig
        end
    end

end


figure;
title('Correct Choice');
histogram(choice_check);
xticks([1,2,3]);
figure;
histogram(rove_check(:,2:3));
hold on;
histogram(rove_check(:,1));
legend('Distractor','Target');
title('Roving');
xticks([0,1]);

