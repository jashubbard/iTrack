function [inst_amp, inst_phase, cum_phase]= hilbert_amp_phase(data)
%This function attempts to convert the result of hilbert output to
%instantaneous phase and amplitudes


%first convert the data with hilbert function
temp_H = hilbert(data);

%then extract real and imaginary part
real_H = real(temp_H);
img_H = imag(temp_H);

%calc instantaneous power!!!
inst_amp = sqrt(power(real_H,2)+power(img_H,2));

%phase time!!
inst_phase = atan2(real_H, img_H);

cum_phase = inst_phase;
for i = 2:1:length(data)
    if inst_phase(i-1) < 0 && inst_phase(i)>0
        cum_phase(i:length(inst_phase)) = cum_phase(i:length(inst_phase))+pi;
    end;
end;

