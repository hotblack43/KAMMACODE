maxshift=0.03
     ig=get_data('ig.txt')
     bg=get_data('bg.txt')
     jd_bg=bg(0,*)
     nbg=N_elements(jd_bg)
     openw,31,'selected_bg.txt'
     for i=0,nbg-1,1 do begin
     idx=where(ig(0,*) ge jd_bg(i)-maxshift and ig(0,*) le jd_bg(i))
     jdx=where(ig(0,*) le jd_bg(i)+maxshift and ig(0,*) ge jd_bg(i))
	if (idx(0) ne -1 and jdx(0) ne -1) then begin
	if (n_elements(idx) ge 1 and n_elements(jdx) ge 1) then printf,31,format='(f15.7,1x,f9.5)',bg(*,i)
	endif
     endfor
     close,31
openw,33,'R.dat'
epsilon=maxshift/1000.
r_best=-1e22
 for shif=-maxshift+epsilon,maxshift-epsilon,0.001 do begin
     ig=get_data('ig.txt')
     bg=get_data('selected_bg.txt')
     ig(0,*)=ig(0,*)+shif
     jd1=max([min(ig(0,*)),min(bg(0,*))])
     jd2=min([max(ig(0,*)),max(bg(0,*))])
     idx=where(ig(0,*) gt jd1 and ig(0,*) lt jd2)
     jd_ig=reform(ig(0,idx));-jd1
     ig=reform(ig(1,idx))
     idx=where(bg(0,*) gt jd1 and bg(0,*) lt jd2)
     jd_bg=reform(bg(0,idx));-jd1
     bg=reform(bg(1,idx))

;    kdx=where(jd_bg lt 6.1 or jd_bg gt 7)
;    jd_bg=jd_bg(kdx)
;    bg=bg(kdx)

     !P.MULTI=[0,1,2]
     plot,jd_ig,ig,xtitle='Days',ytitle='BG (red) and IG'
     oplot,jd_bg,bg,psym=7,color=fsc_color('red')
     xyouts,/normal,0.2,0.9,'Shift on ig time axis: '+string(shif*24.*60,format='(f5.1)')+' min.'
     ig_int=interpol(ig,jd_ig,jd_bg)
     plot,xrange=[0,20],yrange=[0,20],bg,ig_int,psym=7,xtitle='BG',ytitle='Interpolated IG'
; considerthe scatter around the y=x line
     delta=bg-ig_int
     sd=stddev(delta)
     print,'SD of scatter around y=x line is: ',SD
     oplot,[0,20],[0,20],linestyle=2
     oplot,[0,20],[-SD*2,20-SD*2],linestyle=1
     oplot,[0,20],[SD*2,20+SD*2],linestyle=1
     gdx=where(abs(delta) gt SD*2)
     print,'Found ',n_elements(gdx),' outliers'
;..........................................
     r=correlate(bg,ig_int)
     if (abs(r) gt r_best) then begin
	openw,55,strcompress('best_R_lag_series.dat')
	openw,56,strcompress('best_R_lag_metadat.dat')
        print,'Writing out ','best_R_lag_series.dat and best_R_lag_metadat.dat'
	printf,56,shif,r 
	for klm=0,n_elements(bg)-1,1 do printf,55,format='(f15.7,3(1x,f9.4))',jd_bg(klm),bg(klm),ig_int(klm),bg(klm)-ig_int(klm)
        close,55
        close,56
        r_best=abs(r)
     endif
     rmse=sqrt(total(delta/bg)^2)/float(n_elements(bg))
     print,shif,'R(IGi,BG)=',R,rmse
     printf,33,shif*24.*60,R
     xyouts,2,17,'R(IGi,BG)='+string(R,format='(f6.3)')
     endfor
 close,33
 data=get_data('R.dat')
 plot,data(0,*),data(1,*),xtitle='Shift [min]',ytitle='R',ystyle=3
; plot outliers
data=get_data('best_R_lag_series.dat')
plot,xstyle=3,ystyle=3,(data(0,*) mod 1)*24,data(3,*),psym=7,xtitle='Hpurs since noon [GMT?]',ytitle='Residual bg-ig_int at best R'
oplot,[!X.crange],[0,0],linestyle=1
 end
