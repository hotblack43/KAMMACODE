spawn," sed 's/\// /g' petersinput.txt | sed 's/\,/ /g' | sed 's/:/ /g' > petersordnedeinput.dat "
data=get_data('petersordnedeinput.dat')

dd=reform(data(0,*))
mm=reform(data(1,*))
yy=reform(data(2,*))
hh=reform(data(3,*))
mi=reform(data(4,*))
value=reform(data(5,*))
jd=julday(mm,dd,yy,hh,mi)
openw,33,'new.txt'
n=n_elements(jd)
print,'n=',n
for i=0,n-1,1 do begin
printf,33,format='(f15.7,1x,f9.4)',jd(i),value(i)
endfor
close,33
print,'Data now formatted as JD,data in file new.txt'
	end

