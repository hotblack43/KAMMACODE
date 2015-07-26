;
data=get_data('large_meals.dat')
x=reform(data(0,*))	; hours since large morning carbs
idx=sort(x)
data=data(*,idx)
x=reform(data(0,*))	; hours since large morning carbs
y=reform(data(1,*))	; ig
z=reform(data(2,*))	; ig - ig_0
!P.MULTI=[0,1,2]
plot,x,y,psym=3,xtitle='Hours since carbs',ytitle='ig',title='Large morning meals'
oplot,x,smooth(y,31,/edge_truncate)
plot,x,z,psym=3,xtitle='Hours since carbs',ytitle='!7D!3ig'
oplot,x,smooth(z,31,/edge_truncate)
end
