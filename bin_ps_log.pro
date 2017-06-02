
; Bins a power spectrum

;********************************************************************************************************

pro bin_ps_log, f_min=f_min, f_max=f_max, freq=freq, ps=ps, stack_err=stack_err, $
		freq_vector_log, ps_vector, pse_vector

  x=freq

  ; Number of bins
   np = n_elements(freq)
   nlogbin=26.d0

; [i] Low freqs: log binning

    logfmin= alog(f_min) ; max freq for 1/f fitting (log)
    logfmax  = alog(f_max)  ; min freq for 1/f fitting (log)
    d_f = (logfmax-logfmin)/nlogbin

    freq_vector_log =  logfmin + d_f/2.0 + findgen(nlogbin)*d_f ; (log)

    ps_vector    = fltarr(nlogbin)
    pse_vector   = fltarr(nlogbin)

    for i=0, nlogbin-1 do begin
       
       cut = where( alog(x) gt (freq_vector_log[i]-d_f/2.0) and $
                    alog(x) le (freq_vector_log[i]+d_f/2.0), ncut )

       if (ncut ge 3) then begin
         ps_vector[i]  = mean( ps[cut] )  ; binned PS
	  ; error in each bin:
	 if keyword_set(stack_err) then begin
	    pse_vector[i] = sqrt(total(stack_err[cut]^2.d0))/ncut
	 endif else begin
	    pse_vector[i] = stddev(ps[cut])/sqrt(ncut)
; 	    pse_vector[i] = stddev(ps[cut])
	 endelse

       endif

    endfor

cut=where(ps_vector ne 0.d0)
freq_vector_log=freq_vector_log[cut]
ps_vector=ps_vector[cut]
pse_vector=pse_vector[cut]

end

;*****************************************************************************************************