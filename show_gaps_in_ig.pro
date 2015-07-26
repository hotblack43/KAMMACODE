delta_lim=2./24.
ig=get_data('ig.txt')
delta=ig(0,*)-shift(ig(0,*),1)
plot,ig(0,*)-shift(ig(0,*),1),ystyle=3,yrange=[-1,2]
n=n_elements(delta)
for i=1,n-1,1 do begin
if (delta(i) gt delta_lim) then begin
caldat,ig(0,i-1),mm,dd,yy,hh,mi,se
print,'large jump',delta(i)*24.,' hours at: ',mm,dd,yy,hh,mi
endif
endfor
end
