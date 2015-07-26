!P.charsize=3
eps=1.0
;!P.MULTI=[0,1,1]
!P.MULTI=[0,3,3]
ig=get_data('ig.txt')
for ishif=1,9,1 do begin
print,'Shift is: ',ishif*5.0,' minutes'
d_ig=ig(1,*)-shift(ig(1,*),ishif)
dt=ig(0,*)-shift(ig(0,*),ishif)
idx=where(abs(dt*24.*60.-ishif*5.) lt eps)
;plot,yrange=[-10,10],xstyle=3,ystyle=3,(ig(0,idx))*24,d_ig(idx),psym=7,title='Shift is: '+string(fix(ishif*5.0))+' minutes',xtitle='Timer',ytitle='Hop i IG'
plot,yrange=[-10,10],xstyle=3,ystyle=3,(ig(0,idx) mod 1)*24,d_ig(idx),psym=3,title='Shift is: '+string(fix(ishif*5.0))+' minutes',xtitle='Timer siden 12.00',ytitle='Hop i IG'
endfor
end
