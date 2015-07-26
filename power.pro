data=get_data('ig.txt')
x=reform(data(0,*))
y=reform(data(1,*))
;y=randomn(seed,n_elements(y))
;y=y-mean(y)
period=1.0
;y=y+10.*sin(x*!pi*2./period)
scale_factor=1.0
jmax=1.0
value=1.0
Result = LNP_TEST( X, Y , /DOUBLE, HIFAC=scale_factor, JMAX=variable, OFAC=value, WK1=f, WK2=powr)
print,result
plot_oo,1/f,powr,xtitle='Period [days]',psym=10,nsum=3
end
