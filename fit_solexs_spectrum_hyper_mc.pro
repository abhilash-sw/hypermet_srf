spec_in_inds = [1,2,5,10,15,25,35,50,65,75]
n_inds = n_elements(spec_in_inds)

;for iii=0,999 do begin
;i_str = string(iii,format='(I04)')
spectrum_dir = '/media/abhilash/bea6a276-1685-41b3-accd-b394aa7c0539/lalitha/suryahome/suryahome/abhilash/chianti_spectrum/solexs_spectrum/hypermet_60eV/coronal/SDD1/';+i_str

rsp_dir = '/media/abhilash/bea6a276-1685-41b3-accd-b394aa7c0539/lalitha/suryahome/suryahome/abhilash/caldbs/hypermet_caldbs/rsps'

for i_iter=0,999 do begin
i_str = string(i_iter,format='(I04)')
openw,1, '/home/abhilash/abhilash/SOLEXS/srf/hypermet_srf/fit_results/vth_fit_results_'+i_str+'.txt'

o = ospex(/no_gui)
o -> set, spex_file_reader='solexs_read'
set_logenv, 'OSPEX_NOINTERACTIVE', '1'

;spawn, 'rm -rf '+ spectrum_dir+'/fit_results'
;spawn, 'mkdir '+ spectrum_dir+'/fit_results'

rsp_file = rsp_dir+'solexs_hypermet_resp_'+i_str+'.rsp'

for i=0,n_inds-1 do begin       ; start temperature loop
    print, i
    jj = spec_in_inds(i)
    iso_logt=6.0+jj*0.02
    pha_file = spectrum_dir + '/coronal_ch10_3e10_1e27_'+string(iso_logt,format='(f4.2)')+'.ascii_SDD1.pha'
    print, pha_file
    o->set, spex_specfile=pha_file
    o->set, spex_drmfile=rsp_file
    o->set, fit_function='vth'
    o->set, fit_comp_param=[1., 10^iso_logt*0.08/1e6, 1.]
    o->set, fit_comp_free = [1,1,0]
    o->set, fit_comp_model=['chianti', '']
    o->set, fit_comp_minima = [1e-20,0.1,0.1]
    o->set, spex_fit_manual=0, spex_autoplot_enable=0, spex_fitcomp_plot_resid=0, spex_fit_progbar=0

    d=o->getdata(class='spex_fit')
    s = o->get(/fit_comp_params)
    sigmas = o -> get(/fit_comp_sigmas)

    ;spawn, 'basename '+ pha_file + ' .pha', basename
    ;o->savefit, outfile= spectrum_dir + '/fit_results/ospex_result_'+ basename(0) + '.fits'

    ;printf,1, iso_logt, s(0), s(1), s(2), sigmas(0), sigmas(1), sigmas(2)
    out_str = string(iso_logt,s,/PRINT)
    out_str = out_str+string(sigmas,/PRINT)
    printf,1, out_str;string(iso_logt,s,sigmas)
endfor

close, 1

;endfor ;iii
endfor
end
