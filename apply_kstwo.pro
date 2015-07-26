; Test of whether two sequences of numbers are likely drawn from the same parent distribution
; aka Kolmogorov-Sminroff test
;
; read in two arrays
;
dat_head=get_data('kammas_infu_insulin_halfdays_head.dat')
dat_tail=get_data('kammas_infu_insulin_halfdays_tail.dat')
; apply KS2 test
kstwo,dat_head,dat_tail,d,prob
print,'D,prob: ',d,prob
print,'prob near 1 means the two samples are very similar'
end
