
; March'14

; For a given observation, reads power spectra for the MFI (32 channels + polarization)
; Fits each powerspectrum with "fit_ps_logps2.pro"
; Output: text file + plot for each horn
;**************************************************************************************

@fit_ps_logps3

; @fit_ps_logps_fixedalpha

pro final1f, folder, savepng=savepng

fmodelo='fmodelo0'
outname=folder+'_'+fmodelo
; outname=folder

; QUIJOTE Info
info_quijote,mfi=mfi

doplot=1
dopng=0
errors_mc=0
silent=0

datafile=outname+'_fit.txt'
openw, 23, datafile

; resolve_routine, 'fit_ps_logps2', /COMPILE_FULL_FILE
; resolve_routine, 'fit_ps_logps3', /COMPILE_FULL_FILE

;**********************
; Single channels
;**********************

FOR horn=1,4 DO BEGIN

;  Select channels for each horn

   ihorn=horn-1
   ch_st=ihorn*8
   ch_end=ihorn*8+7


  window, 19, xsize=800, ysize=900, tit=folder
  !p.multi=[0,2,4]

  for j=ch_st,ch_end do begin

  chan=j+1

  print, " "
  print," >> Working on channel", chan

          fit_ps_logps3, folder=folder, chan=chan, $
	  doplot=doplot, dopng=dopng, errors_mc=errors_mc, silent=silent, $
	  fmodelo=fmodelo, /peakrm, a_in, a_out, sigma, redchisq, c, sigma_mc
 
	  printf, 23, mfi.chan_name[chan-1], c, redchisq, a_out[0], sigma[0], $
							 a_out[1], sigma[1], $
; 							 a_out[2], sigma[2], $
; 							 a_out[3], sigma[3], $
							 FORMAT='(A5, I3, F6.1, 4F8.2, 4F7.3)'

  endfor
 !p.multi=0

	if keyword_set(savepng) then begin

	horn_id=strtrim(horn,2)

    	outfile=outname+'_horn'+horn_id+'.png'
    	im = tvrd(true=1)
    	write_png, outfile, im, $
    	xresolution=1800, yresolution=1200


	out_dir='/net/trevina/scratch/alba/noise1f/'+folder+'/images/'
	test1=file_test(out_dir, /directory)
	if not (test1) then file_mkdir, out_dir
	file_move, outfile, out_dir, /overwrite

	endif

ENDFOR

;**********************
; Polarization channels
;**********************

chanid=['111','113','217','219','311','313','417','419']

window, 17, xsize=800, ysize=900, tit=folder
!p.multi=[0,2,4]

for j=0,7 do begin

  chan=chanid[j]

  print, " "
  print," >> Working on polarization channel ", chanid[j]

          fit_ps_logps3, folder=folder, chan=chan, $
	  doplot=doplot, dopng=dopng, errors_mc=errors_mc, silent=silent, $
	  fmodelo=fmodelo, /peakrm, a_in, a_out, sigma, redchisq, c, sigma_mc

	  printf, 23, chan, c, redchisq, a_out[0], sigma[0], $
			       a_out[1], sigma[1], $
			       a_out[2], sigma[2], $
; 			       a_out[3], sigma[3], $
			       FORMAT='(A5, I3, F6.1, 4F8.2, 4F7.3)'

endfor
!p.multi=0

	close, 23
        out_dir='/net/trevina/scratch/alba/noise1f/'+folder
	file_move, datafile, out_dir, /overwrite

	if keyword_set(savepng) then begin

	horn_id=strtrim(horn,2)

    	outfile=outname+'_pol.png'
    	im = tvrd(true=1)
    	write_png, outfile, im, $
    	xresolution=1800, yresolution=1200

	out_dir='/net/trevina/scratch/alba/noise1f/'+folder+'/images/'
	file_move, outfile, out_dir, /overwrite

	endif

end
