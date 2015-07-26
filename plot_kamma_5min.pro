PRO ploit,data,ix,iy
names=['timestamp','ig','bg','ba','bo','ca','is','ss','ill']
plot,xstyle=3,ystyle=3,data(ix,*),data(iy,*),xtitle=names(ix),ytitle=names(iy),charsize=1.7,psym=7
return
end
data=get_data('./kamma_wrap/julian-dates/kamma_5_minute.dat')
idx=where(data(1,*) ne -1)
help,data
data=data(*,idx)
help,data
;
jd=reform(data(0,*))
ig=reform(data(1,*))
bg=reform(data(2,*))
ba=reform(data(3,*))
bo=reform(data(4,*))
ca=reform(data(5,*))
is=reform(data(6,*))
ss=reform(data(7,*))
ill=reform(data(8,*))
;
ploit,data,1,8

end















