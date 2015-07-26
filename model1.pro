PRO insulin_range,ca,g0,frd0,x0min,x0max
  EPS=1.7 ; upper bound on the size of the deviation to target T after 4H
  ; sample regression model
  const= 46.111146
;       4.1955383 ±        38.092438   ; ca
;      -5639.7037 ±        2022.6390   ; babo
;      -40.157351 ±        23.886793   ; frd
;     -0.58942678 ±       0.21538708   ; ig0
  coef=[ -5639.7037, 4.1955383, -40.157351, -0.58942678 ]
  x0min = 1/coef(0)*((g0-eps)-const-coef(1)*ca-coef(2)*frd0-g0*coef(3))
  x0max = 1/coef(0)*((g0+eps)-const-coef(1)*ca-coef(2)*frd0-g0*coef(3))
end

g0init=4.5
frd0=0.75                   ; breakfast time 
for ca=15,25,5 do begin     ; loop over carb sizes
  for mi=0,300,30 do begin  ; loop over meal init time
    frd=frd0+mi/60./24.     ; julian date
    for i=0,6 do begin      ; loop over ig0
      g0=g0init+0.5*i
      insulin_range,ca,g0,frd,x0min,x0max
      print,frd,g0,ca,x0min,x0max,mean([ca/(x0max*12.),ca/(x0min*12)]) ; IC-setting
    endfor
  endfor
endfor
end
