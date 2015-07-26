path='./kamma_wrap/julian-dates/'
ig=get_data(path+'ig.txt')
bg=get_data(path+'bg.txt')
ca=get_data(path+'carbs.txt')
bo=get_data(path+'bolus.txt')
ba=get_data(path+'basal.txt')
is=get_data(path+'infusion_shift.txt')      
ss=get_data(path+'sensorshift.txt')      
ill=get_data(path+'illness.txt')      
openw,44,'kamma_5_minute.dat'
;for ijd=min(long(ig(0,*))),min(long(ig(0,*)))+2,1 do begin
for ijd=min(long(ig(0,*))),max(long(ig(0,*)))-1,1 do begin
  print,'Day: ',ijd
  caldat,ijd,mm,dd,yy
  idx_ill=where(long(ill(0,*)) eq ijd)
  if (idx_ill(0) ne -1) then begin
    p_ill=1.0
  endif else begin
    p_ill=0.0
  endelse
  for hour=0,23,1 do begin
    for minu=0,55,5 do begin
      p_ig=-1.0
      p_bg=-1.0
      p_ca=0.0
      p_bo=0.0
      p_ba=0.0
      p_is=0.0
      p_ss=0.0
      ijdts=julday(mm,dd,yy,hour+12,minu, 0)
      ijdtsn=julday(mm,dd,yy,hour+12,minu+4.9, 0)
;      print,format='(3(f15.7,1x))', ijd,ijdts,ijdtsn
      idx_ig=where(ig(0,*) ge ijdts and ig(0,*) lt ijdtsn)
      idx_ca=where(ca(0,*) ge ijdts and ca(0,*) lt ijdtsn)
      idx_ba=where(ba(0,*) ge ijdts and ba(0,*) lt ijdtsn)
      idx_bo=where(bo(0,*) ge ijdts and bo(0,*) lt ijdtsn)
      idx_bg=where(bg(0,*) ge ijdts and bg(0,*) lt ijdtsn)

      ; find relevant infusionset ID
      is_tdiff=ijdts-is     
      idx_is=where(is_tdiff ge 0.0)
      last_is=min(is_tdiff(idx_is)) 
      idx_isfound=where(is(0,*) le ijdts-last_is)
      if (idx_isfound(0) ne -1) then begin
        p_is=max(is(1,idx_isfound-1))
      endif

      ; find relevant sensor ID
      ss_tdiff=ijdts-ss     
      idx_ss=where(ss_tdiff ge 0.0)
      last_ss=min(ss_tdiff(idx_ss)) 
      idx_ssfound=where(ss(0,*) le ijdts-last_ss)
      if (idx_ssfound(0) ne -1) then begin
        p_ss=max(is(1,idx_ssfound-1))
      endif

      if (idx_bg(0) ne -1) then begin
        p_bg=bg(1,idx_bg)
        if (n_elements(idx_bg) gt 1) then begin
          print,'WARNING - multiple BG values - using mean'
          caldat,ijdts,mm,dd,yy,hh,mi
          print,ijdts,mm,dd,yy,hh,mi
          p_bg=mean(bg(1,idx_bg))
        endif
      endif

      if (idx_ig(0) ne -1) then begin
        p_ig=ig(1,idx_ig)
        if (n_elements(idx_ig) gt 1) then begin
;          stop
          caldat,ijdts,mm,dd,yy,hh,mi
          print,'WARNING - multiple IG values - should not be possible - skip them'
          print,ijdts,mm,dd,yy,hh,mi
          p_ig=-1.0
        endif
      endif
      if (idx_ca(0) ne -1) then begin
        p_ca=total(ca(1,idx_ca))
      endif
      if (idx_ba(0) ne -1) then begin
        p_ba=total(ba(1,idx_ba))
      endif
      if (idx_bo(0) ne -1) then begin
        p_bo=total(bo(1,idx_bo))
      endif
      printf,44,format='(1(f15.7,1x),8(f7.2,1x))',ijdts,p_ig,p_bg,p_ba,p_bo,p_ca,p_is,p_ss,p_ill
    endfor
  endfor
endfor
close,44
end

