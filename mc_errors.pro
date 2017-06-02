

pro mc_errors, fmodelo, freq, ps, ps_err, sigma_mc=sigma_mc, nsimu_fin=nsimu_fin, doplot=doplot

if not keyword_set(doplot) then doplot=0

nsimu=1000.d0

n=n_elements(ps)

if fmodelo eq 'fmodelo4' then nparam=4 else nparam=3

;Storage

est_ps=dblarr(n, nsimu)

mean_ps=dblarr(n)
sigma_ps=dblarr(n)

est_a=dblarr(nparam,nsimu)
est_a_sigma=dblarr(nparam,nsimu)

count=0.d0

; Begin simulations
for i=0, nsimu-1L do begin

ps_sim_err=randomn(seed,n)*ps_err
ps_sim=ps+ps_sim_err

est_ps[*,i]=ps_sim

fit_ps_logps2, freq=freq, powerspec=ps_sim, ps_err=ps_err, doplot=0, dopng=0, $
        errors_mc=0, fmodelo=fmodelo, silent=1, a_in, a_out, sigma, redchisq, c

;   Consider the simulation only in case that the fit has converged

    if c eq 1L then begin

    est_a[*,i]=a_out
    est_a_sigma[*,i]=sigma

    endif else begin

    count=count+1.d0

    endelse

endfor

nsimu_fin=nsimu-count
cut=where(est_a[0,*] ne 0.0)

; Mean parameters

mean_a=total(est_a[0,cut])/nsimu_fin
stddev_a=sqrt(total((est_a[0,cut]-mean_a)^2.)/nsimu_fin)

mean_fk=total(est_a[1,cut])/nsimu_fin
stddev_fk=sqrt(total((est_a[1,cut]-mean_fk)^2.)/nsimu_fin)

mean_alpha=total(est_a[2,cut])/nsimu_fin
stddev_alpha=sqrt(total((est_a[2,cut]-mean_alpha)^2.)/nsimu_fin)

if nparam eq 4 then begin
mean_beta=total(est_a[3,cut])/nsimu_fin
stddev_beta=sqrt(total((est_a[3,cut]-mean_beta)^2.)/nsimu_fin)
endif else begin
mean_beta=''
stddev_beta=''
endelse

; Mean power spectrum and error

for j=0, n-1 do begin
mean_ps[j]=total(est_ps[j,cut])/nsimu_fin
sigma_ps[j]=sqrt(total((est_ps[j,cut]-mean_ps[j])^2.)/nsimu_fin)
endfor

sigma_mc=[stddev_a, stddev_fk, stddev_alpha, stddev_beta]

; Print output

	print, "****************************************************************************"
	print, "  > MC: Output parameters - ", mean_a, mean_fk, mean_alpha, mean_beta
	print, "  > MC: Errors            - ", stddev_a, stddev_fk, stddev_alpha, stddev_beta
	print, "  > Total number of valid simulations - ", uint(nsimu_fin)
	print, "****************************************************************************"
	print, "****************************************************************************"

; Plots

if doplot then begin

    window, 2
    kk1=sigma_ps/mean_ps
    ; plot, b.freq, kk1, /xs,/ys,/xlog,xrange=[0.01,500.d0], yrange=[0.01,0.6]
    plot_oo, b.freq, mean_ps, /xs,/ys,/xlog,xrange=[0.01,500.d0],yrange=[1.d-9,1.d-2]         
    errplot, b.freq, mean_ps-sigma_ps, mean_ps+sigma_ps

    window,3
    kk2=b.sigma/b.ps
    ; plot, b.freq, kk2, col=fsc_color('red'),/xs,/ys,/xlog,xrange=[0.01,500.d0], yrange=[0.01,0.6]
    plot_oo, b.freq, b.ps, col=fsc_color('red'),/xs,/ys,/xlog,xrange=[0.01,500.d0],yrange=[1.d-9,1.d-2]   
    errplot,b.freq,b.ps-b.sigma,b.ps+b.sigma, col=fsc_color('red')

    window,4

    kk3=sigma_ps/b.sigma
    plot, b.freq, kk3, col=fsc_color('red'),/xs,/ys,/xlog,xrange=[0.01,500.d0], yrange=[0.01,6]
    ; oplot, freq_vector, ps_vector, col=fsc_color('red')         
    ; errplot,freq_vector, ps_vector-pse_vector, ps_vector+pse_vector, col=fsc_color('red')

endif

end