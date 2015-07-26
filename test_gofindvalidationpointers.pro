PRO gofindvalidationpointers,nrows,idx,jdx
; find the integers not present in the list idx between 0 and nrows-1, return in jdx
; allow repeats in jdx as in idx
list=[]
for k=0,nrows-1,1 do begin
kdx=where(idx eq k)
if (kdx(0) eq  -1) then list=[list,k]
endfor
; 'list' is now a list of the integers NOT present in idx, now bootstrap them, with replacement
ldx=long(randomu(seed,n_elements(idx))*n_elements(idx))
jdx=list(ldx)
print,idx
print,jdx
stop
return
end

nrows=12
bootfraction=0.5
idx=long(randomu(seed,nrows*bootfraction)*nrows)
gofindvalidationpointers,nrows,idx,jdx
end
