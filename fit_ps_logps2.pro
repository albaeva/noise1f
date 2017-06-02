;
; 9/8/2013
; Modified: February 2014

; Fits a noise power spectrum saved as a fits file. Correlated or normal channels.

; obs_name = 

;	"synth" --> synthetic pink noise generated with "noisegen.pro"
;	"COSMO2ELLIP-130223-2315" ; fixed telescope; elliptic filter on
;	"COSMO2ELLIP-130226-0120" ; tracking on; elliptic filter on
;	"COSMO2A-130313-2319" ; fixed telescope; elliptic filter on (long observation)
;	"MOONFON121206_0429" ;telescope raster; elliptic filter on
;	"COSMO3D-130729-0002"; 2.5 h observation
;	"COSMO31FB-131001-1435" ; fixed telescope (axis not disabled)


; dopng = > save png (0:no ; 1:yes)

;----------------------------------------------------------------------------

@mc_errors
@bin_ps_log

;---------------------------------------------------------------------------------
function fmodelo1, y, a

x=exp(y)

; ;Model function for noise PS:

fun=a[0]*(1.d0+(a[1]/x)^a[2])
log_fun=alog(fun)
; 
;Derivatives with respect to the parameters:

df0=(1.d0+(a[1]/x)^a[2])/fun
df1=((a[0]*a[2]/a[1])*(a[1]/x)^a[2])/fun
df2=a[0]*((a[1]/x)^a[2])*alog(a[1]/x)/fun
 
return, [[log_fun],[df0],[df1],[df2]]
end 
;---------------------------------------------------------------------------------
function fmodelo2, y, a

x=exp(y)

; ;Model function for noise PS:
; 
fun=a[0]*(1.d0+(a[1]/x)^a[2])^2.d0
log_fun=alog(fun)

; ;Derivatives with respect to the parameters:
; 
df0=((1.d0+(a[1]/x)^a[2])^2.d0)/fun
df1=(2.d0*(a[0]*a[2]/a[1])*(1.d0+(a[1]/x)^a[2])*(a[1]/x)^a[2])/fun
df2=(2.d0*a[0]*(1.d0+(a[1]/x)^a[2])*((a[1]/x)^a[2])*alog(a[1]/x))/fun
 
return, [[log_fun],[df0],[df1],[df2]]
end 
;---------------------------------------------------------------------------------
function fmodelo3, y, a

x=exp(y)

; ;Model function for noise PS:
; 
fun=a[0]*(1.d0+(a[1]/x))^a[2]
log_fun=alog(fun)
; 
; ;Derivatives with respect to the parameters:
; 
df0=((1.d0+(a[1]/x))^a[2])/fun
df1=(a[0]*a[2]*((1.d0+(a[1]/x))^(a[2]-1.d0))/x)/fun
df2=(a[0]*((1.d0+(a[1]/x))^a[2])*alog(1.d0+(a[1]/x)))/fun
 
return, [[log_fun],[df0],[df1],[df2]]
end 
;---------------------------------------------------------------------------------
function fmodelo4, y, a

x=exp(y)

; ;Model function for noise PS:

rind=a[2]+a[3]*y
fun=a[0]*(1.d0+(a[1]/x)^rind)
log_fun=alog(fun)
; 
;Derivatives with respect to the parameters:

df0=(1.d0+(a[1]/x)^rind)/fun
df1=((a[0]*rind/a[1])*(a[1]/x)^rind)/fun
df2=a[0]*((a[1]/x)^rind)*alog(a[1]/x)/fun
df3=a[0]*((a[1]/x)^rind)*alog(a[1]/x)*y/fun
 
return, [[log_fun],[df0],[df1],[df2],[df3]]

end 
;---------------------------------------------------------------------------------


pro fit_ps_logps2, freq=freq, powerspec=powerspec, ps_err=ps_err, $
    folder=folder, ff_fits=ff_fits, chan=chan, $
    doplot=doplot, dopng=dopng, errors_mc=errors_mc, silent=silent, $
    fmodelo=fmodelo, peakrm=peakrm, a_in, a_out, sigma, redchisq, c, sigma_mc

color8

if not keyword_set(doplot) then doplot=0
if not keyword_set(dopng) then dopng=0
if not keyword_set(errors_mc) then errors_mc=0
if not keyword_set(silent) then silent=0

; QUIJOTE Info

	info_quijote,mfi=mfi
	fsamp=1.d3

; [1] If no input power spectrum, read file

 if not keyword_set(folder) and not keyword_set(powerspec) and not keyword_set(ff_fits) then begin
    print, 'Power spectra folders'
    print, '---------------------'
    spawn, 'ls /net/trevina/scratch/alba/noise1f/* -d | xargs -n1 basename' 
    folder=''
    read, 'Enter folder name: ', folder
    chan=0
    read, 'Enter channel number (1-32, 111, 113, 217, 219, 311, 313, 417, 419): ', chan
 endif

 if keyword_set(folder) then begin
    if chan le 32 then chan_id=mfi.chan_name[chan-1] else chan_id=strtrim(string(chan),2)
    path2fits = "/net/trevina/scratch/alba/noise1f/"+folder+'/'
;   Observation name
    ind_bin=stregex(folder, '-v')
    if ind_bin ne -1 then obs_name=strmid(folder,0,ind_bin) else obs_name=folder
    ff_fits=path2fits+obs_name+'_chan'+chan_id+'.fits' 
 endif

 if keyword_set(ff_fits) then begin
        a=mrdfits(ff_fits, 1, hdr) 
        freq=a.freq
	powerspec=a.ps*1.d12 ; in microV^2*s
	if tag_exist(a, 'stack_err') then stack_err=a.stack_err*1.d12 ; in microV^2*s
	tbin = FXPAR( hdr, 'TBIN(S)', COUNT=count1) 
	tbase = FXPAR( hdr, 'TBASE(S)', COUNT=count2)
 endif else begin
 count1=0 & count2=0
 endelse 

	n=n_elements(freq)
        fcut = freq[n-1]
 
; [2] Range of frequencies to fit

  	if count1 eq 1 then f_max=fcut else f_max=350.d0
	if count2 eq 1 then begin
	    f_min=1./tbase 
	endif else begin
	    f_min=0.001  ;1.d0/(n/1000.)
	endelse

; [3] Remove high peak data

      if keyword_set(peakrm) then peak_rm, freq, powerspec, f_max, stack_err, silent=silent

; [5] If input power spectrum is not yet binned, bin it

      if not keyword_set(ps_err) then begin

	  if not silent then print,"  > Binning power spectrum"

; 	  Bin power spectrum
	  if keyword_set(stack_err) then begin
	  bin_ps_log, f_min=f_min, f_max=f_max, freq=freq, ps=powerspec, stack_err=stack_err, $
		freq_vector_log, ps_vector, pse_vector
	  endif else begin
	  bin_ps_log, f_min=f_min, f_max=f_max, freq=freq, ps=powerspec, $
		freq_vector_log, ps_vector, pse_vector
	  endelse
 
;    [6] Convert PS to log

	  ps_vector_log=alog(ps_vector)
	  pse_vector_log=pse_vector/abs(ps_vector)

	endif else begin

	  freq_vector_log=freq
	  ps_vector_log=powerspec
	  pse_vector_log=ps_err

      endelse

; [7] Estimate parameters

  estimate_param, fmodelo, freq_vector_log, ps_vector_log, pse_vector_log, a_in, silent=silent
  
; [8] Fit full power spectrum

  full_fit, freq_vector_log, ps_vector_log, pse_vector_log, $
            a_in, a_out, sigma, redchisq, c, fmodelo=fmodelo, silent=silent

; [9] Determine errors through Montecarlo simulations

  if errors_mc then begin

   if not silent then print,"  > Calculating MC errors..."

   mc_errors, fmodelo, freq_vector_log, ps_vector_log, pse_vector_log, sigma_mc=sigma_mc, $
	      nsimu_fin=nsimu_fin, doplot=0

  endif

; [10] Sensitivity

; sig_v=sqrt(a_out[0])
; err_sig_v=sigma[0]/(2.d0*sig_v)

; [11] Plot

if doplot then begin

   x=exp(freq_vector_log)

   ymin=min(powerspec[where(freq lt f_max)])*100000.d0
   ymax=max(powerspec[where(freq gt f_min)])

;   ymin=min(ps_vector[where(freq_vector_log lt f_max)])
;   ymax=max(ps_vector[where(freq_vector_log gt f_min)])

;      ymax=1.d8
;      ymin=10.d0

;   if keyword_set(folder) then plottit= obs_name + ' - channel '+chan_id else plottit=''
   if keyword_set(folder) then plottit= 'channel '+chan_id else plottit=''

  plot, freq, powerspec, /xs,/ys,/xlog,/ylog,xrange=[0.002,fcut], yrange=[ymin, ymax], $
  tit= plottit, xtit='frequency (Hz)',ytit='PS (!4l!XV!u2!ns)', $ 
  charsize=2.d0, background=255, color=0

  oplot, x, ps_vector, psym=2,col=2
  errplot, x, ps_vector-pse_vector, ps_vector+pse_vector,col=2

  case fmodelo of
 'fmodelo1' : oplot, x, exp(fmodelo1(freq_vector_log, a_in)) , col=4, thick=1.5
 'fmodelo2' : oplot, x, exp(fmodelo2(freq_vector_log, a_in)) , col=4, thick=1.5
 'fmodelo3' : oplot, x, exp(fmodelo3(freq_vector_log, a_in)) , col=4, thick=1.5
 'fmodelo4' : oplot, x, exp(fmodelo4(freq_vector_log, a_in)) , col=4, thick=1.5
  endcase

  case fmodelo of
 'fmodelo1' : oplot, x, exp(fmodelo1(freq_vector_log, a_out)) , col=3, thick=2.0
 'fmodelo2' : oplot, x, exp(fmodelo2(freq_vector_log, a_out)) , col=3, thick=2.0
 'fmodelo3' : oplot, x, exp(fmodelo3(freq_vector_log, a_out)) , col=3, thick=2.0
 'fmodelo4' : oplot, x, exp(fmodelo4(freq_vector_log, a_out)) , col=3, thick=2.0
  endcase

  a_wn='!4r!X!u2!n='+string(a_out[0], FORMAT='(F6.2)')+'!9+!X'+string(sigma[0], FORMAT='(F5.2)')
  f_k='f!dk!n='+string(a_out[1], FORMAT='(F6.2)')+'!9+!X'+string(sigma[1], FORMAT='(F6.2)')
  alpha='!4a!3='+string(a_out[2], FORMAT='(F6.3)')+'!9+!X'+string(sigma[2], FORMAT='(F6.3)')
  if fmodelo eq 'fmodelo4' then beta='!4b!X='+string(a_out[3], FORMAT='(F6.3)')+$
                           '!9+!X'+string(sigma[3], FORMAT='(F6.3)')
  chi2='!4V!3!u2!n='+string(redchisq, FORMAT='(F5.1)')

;   xyouts, 95, 180, alpha, /device, charsize=1.4, col=0in
;   xyouts, 95, 160, f_k, /device, charsize=1.4, col=0
;   xyouts, 95, 140, chi2, /device, charsize=1.4, col=0

  ytic=(alog(ymax)-alog(ymin))/10.d0

  xyouts, 4.d0, ymin*exp(9.*ytic), a_wn, charsize=1.3, col=0
  xyouts, 4.d0, ymin*exp(8.*ytic), f_k, charsize=1.3, col=0
  xyouts, 4.d0, ymin*exp(7.*ytic), alpha, charsize=1.3, col=0
  if fmodelo eq 'fmodelo4' then xyouts, 4.d0, ymin*exp(6.*ytic) , beta, charsize=1.3, col=0
  if not keyword_set(chan) then chan=-1
  if (chan le 32) then begin
     yprintc=ymin*exp(1.5*ytic) 
     yprintchi=ymin*exp(0.5*ytic)
  endif else begin
     yprintc=ymin*exp(9.*ytic) 
     yprintchi=ymin*exp(8.*ytic)
  endelse
  if (c eq 1) then xyouts, 0.003, yprintc, 'Converged', charsize=1.3, col=0
  xyouts, 0.003, yprintchi, chi2, charsize=1.3, col=0
; 
;   legend, ['True PS', 'Binned PS', 'Estimated fit', 'Full fit'], psym=0 , col=[0,2,4,3], $
;   /bottom, /left, textcolor=0, box=0, charsize=1.2

  if dopng then begin

    ; Save png

    outfile=''

    read, 'Save png image as... : ', outfile

  ;   outfile='chan_'+chan_id+'_fit.png'

    if outfile ne '' then begin

	im = tvrd(true=1)
	write_png, outfile, im, $
	xresolution=1800, yresolution=1200

    ; Move png to directory

;     if keyword_set(folder) then begin
; 
; 	out_dir='/net/trevina/scratch/alba/noise1f/'+folder+'/images/'
; 
; 	test1=file_test(out_dir, /directory)
; 
; 	if not (test1) then file_mkdir, out_dir
;     
; 	file_move, outfile, out_dir, /overwrite
; 
;     endif

    endif

    endif
    
endif


end


;*********************************************************************************************

pro peak_rm, freq, ps, f_max, stack_err, silent=silent

    n=f_max/50.d0
    kk=(lindgen(n)+1L)*50L
    tot=0L
    for i=0L, n-1L do begin
      ind=where(long(freq) eq kk[i], count, complement=good)
      tot=tot+count
      freq=freq[good]
      ps=ps[good]
      if keyword_set(stack_err) then stack_err=stack_err[good]
;       ind=where(long(freq) eq kk[i]-1.d0, count, complement=good)
;       tot=tot+count
;       freq=freq[good]
;       ps=ps[good]
    endfor

  if not silent then print,"  > Removed ", tot, " high peak data"

end

;*********************************************************************************************

pro estimate_param, fmodelo, freq_vector_log, ps_vector, pse_vector, a_in, silent=silent

  ; Cut between high and low binning:

  bincut= 1.d0
  log_bincut=alog(bincut)

  low=where(freq_vector_log le log_bincut, complement=high)

  lf_vector=freq_vector_log[low]
  hf_vector=freq_vector_log[high]
  lps_vector=ps_vector[low]
  hps_vector=ps_vector[high]
  lpse_vector=pse_vector[low]
  hpse_vector=pse_vector[high]

  nbintot=n_elements(freq_vector_log)

  if fmodelo eq 'fmodelo4' then a=fltarr(4)  else a=fltarr(3); storage

;   [i] White noise level (estimated from last bin in linear binning)

      nbin=n_elements(hps_vector)
     
      log_wn = hps_vector[nbin-1] ; variance of gaussian noise
      a[0]=exp(log_wn)

;   [ii] Slope   

;     Linear fit in log binning

      if (n_elements(lps_vector) ge 3) then begin
      coef = poly_fit(lf_vector,lps_vector,1, /double, measure_error=lpse_vector)
      endif else begin
      stop, 'Not enough 1/f points for fit!'
      endelse

      if fmodelo eq 'fmodelo2' then a[2] = - coef[1]/2.d0 else a[2] = - coef[1] ; slope

;   [iii] Knee frequency

      kk1=where(coef[0]+coef[1]*freq_vector_log le log_wn)
      if n_elements(kk1) ge 2 then tiruri=1 else tiruri=0
      a[1]=freq_vector_log[kk1[tiruri]] ; knee frequency

      if fmodelo eq 'fmodelo4' then a[3]=-0.05d0

      if not silent then print, "  > Fit: Input parameters  - ", a

      a_in = reform(a)

end

;***************************************************************************************************

pro full_fit, freq_vector, ps_vector, pse_vector, a_in, a_out, sigma, redchisq, c, $
    fmodelo=fmodelo, silent=silent

    if not keyword_set(silent) then silent=0
    nbin=n_elements(ps_vector)
    errors=fltarr(nbin)+1.d0
    a=a_in
    yfit=lmfit(freq_vector,ps_vector,a,function_name=fmodelo,chisq=chisq, $
    sigma=sigma,iter=iter,convergence=c, $
    measure_errors=errors,itmax=100,/double)
    redchisq = chisq/(nbin - n_elements(a))
    a_out = a
 
    if not silent then begin

	print, "  > Fit: Output parameters - ", a_out
	print, "  > Fit: Errors            - ", sigma
	print, "  > Reduced chi squared    - ", redchisq

	if (c eq 1) then print, "  > Converged! " else print, "  > Didn't converge!"

    endif
end

;***************************************************************************************************
