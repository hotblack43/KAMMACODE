

 frd=findgen(100)/99.*(1-0.75)+0.75
 tsc=findgen(100)/99.*(4./24.-0)+0.0
	z=fltarr(100,100)
	ig0=8.0
	for i=0,99,1 do begin
	for j=0,99,1 do begin
	z(i,j)= 124.82991*tsc(j) - 3.47726*ig0 + 315.08064*frd(i) - 264.48277
	endfor
	endfor
	levs=[-400,-300,-200,-100,-50,-25,0,25,50,100,200]
	nlevs=n_elements(levs)
	cthick=fltarr(nlevs)*0+1
	cthick(where(levs eq 0))=4
	clabels=cthick*0+1
	
	!P.CHarsize=2	
	contour,z,frd,tsc,levels=levs,c_thick=cthick,/downhill,xtitle='frd',ytitle='tsc',xstyle=3,ystyle=3,title='diff(ig,insulin), units=[ig units/Insulin units]',c_labels=clabels
	end

