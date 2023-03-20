#####################################################
# @Author: Abhilash Sarwade
# @Date:   2021-10-20 14:19:39
# @email: sarwade@ursc.gov.in
# @File Name: hypermet_simulation.py
# @Project: hypermet_srf

# @Last Modified time: 2023-03-16 09:38:29
#####################################################

"""
Four Possible Edge Filters
1. Cu, Ni (Target, Filter)
2. Co, Fe
3. Fe, Mn
4. Mo, Zr
"""

import numpy as np
import solexs_caldbgen
from solexs_caldbgen.calc_srf import gaussian, calc_xsm_sigma
from scipy.special import erf, erfc
from scipy import optimize
from xraydb import XrayDB
import corner

xdb = XrayDB()

def calc_xsm_i_esc(E_ph,p,q,r):
     return p*(np.exp(-q*(E_ph-1.739))+r)

def hypermet_func(E_ph_array,ene_out_512,a,b,c,i_tail0,i_shelf,alpha,beta,gamma,p,q,r): #E_ph in keV
    si_ka = 1739 #eV
    enes_out = ene_out_512#np.loadtxt('../../ebounds/ebounds_out/energy_bins_512_out.dat')
    enes_out_max = np.max(enes_out)
    enes_out_mean = (enes_out[:,0]+enes_out[:,1])/2

    n_ch = len(enes_out)#512
    del_e =enes_out[0,1] - enes_out[0,0] #60e-3 #keV (assuming uniform bin widths)
    ch = np.arange(n_ch)

    try:
        n_E = len(E_ph_array)
    except:
        E_ph_array = [E_ph_array]

    spec = np.zeros((len(E_ph_array),n_ch))

    for iii,E_ph in enumerate(E_ph_array):
        # if E_ph>enes_out_max:
        #     # return spec
        #     continue

        spectrum_p1 = np.zeros(n_ch)
        spectrum_p2 = np.zeros(n_ch)
        spectrum_t = np.zeros(n_ch)
        spectrum_s = np.zeros(n_ch)

        ch_ph = int(E_ph/del_e)

        sigma_main = calc_xsm_sigma(E_ph,a,b,c)

        i_tail = i_tail0*(E_ph/15)**gamma

        spectrum_p1 = gaussian(enes_out_mean,E_ph,sigma_main)
        spectrum_t = i_tail*np.exp((enes_out_mean-E_ph)/beta)*erfc((enes_out_mean-E_ph)/np.sqrt(2)/sigma_main + np.sqrt(2)*sigma_main/2/beta)
        spectrum_s = i_shelf*(enes_out_mean/10)**(-1*alpha)*(erf((enes_out_mean - si_ka/1000)/np.sqrt(2)/sigma_main) + 2)
        spectrum_s[enes_out_mean>E_ph]=0

        if E_ph*1e3 > si_ka:
            E_esc = E_ph - si_ka/1000
            ch_esc = int(E_esc/del_e)

            sigma_esc = calc_xsm_sigma(E_esc,a,b,c)
            i_esc = calc_xsm_i_esc(E_ph,p,q,r)
            spectrum_p2 = i_esc*gaussian(enes_out_mean,E_esc,sigma_esc)

        spectrum_full = spectrum_p1 + spectrum_p2 + spectrum_t + spectrum_s
        spec[iii,:] = spectrum_full/np.sum(spectrum_full)

    return spec



##### Simulation
a=3528.9
b=299.04
c=9.41
i_tail0 = 6.21e-3
gamma = 4.63
beta = 0.5
i_shelf = 1e-4
alpha = 0.3

p = 0.017
q = 0.4189
r = 0.04065



esc_peak_int_arr = np.loadtxt('/home/abhilash/abhilash/SOLEXS/caldbs/hypermet_60eV/aditya-l1/solexs/data/bcf/electronic_noise/xsm_escape_peak_intensity.dat')
esc_peak_int_arr[np.isnan(esc_peak_int_arr[:,1]),1] = 0

ene_out_512 = np.loadtxt('/home/abhilash/abhilash/SOLEXS/caldbs/hypermet_60eV/aditya-l1/solexs/data/cpf/ebounds/ebounds_out/energy_bins_512_out.dat')

count_rate = 3 # counts per second
time_int = 24*3*3600 # 3 days
total_counts = count_rate*time_int 

si_ka = 1739 #eV

E_fe = xdb.xray_lines('Fe')['Ka1'].energy/1e3
E_co = xdb.xray_lines('Co')['Ka1'].energy/1e3
E_cu = xdb.xray_lines('Cu')['Ka1'].energy/1e3
E_mo = xdb.xray_lines('Mo')['Ka1'].energy/1e3


#Fe
# sigma_main = calc_xsm_sigma(E_fe)
# i_tail = i_tail0*(E_fe/15)**gamma
# i_esc = np.interp(E_fe,esc_peak_int_arr[:,0],esc_peak_int_arr[:,1])
# E_esc = E_fe - si_ka/1000
# sigma_esc = calc_xsm_sigma(E_esc)

spec_fe = hypermet_func(E_fe,ene_out_512,a,b,c,i_tail0,i_shelf,alpha,beta,gamma,p,q,r)*total_counts

#Co
# sigma_main = calc_xsm_sigma(E_co)
# i_tail = i_tail0*(E_co/15)**gamma
# i_esc = np.interp(E_co,esc_peak_int_arr[:,0],esc_peak_int_arr[:,1])
# E_esc = E_co - si_ka/1000
# sigma_esc = calc_xsm_sigma(E_esc)

spec_co = hypermet_func(E_co,ene_out_512,a,b,c,i_tail0,i_shelf,alpha,beta,gamma,p,q,r)*total_counts

#Cu
# sigma_main = calc_xsm_sigma(E_cu)
# i_tail = i_tail0*(E_cu/15)**gamma
# i_esc = np.interp(E_cu,esc_peak_int_arr[:,0],esc_peak_int_arr[:,1])
# E_esc = E_cu - si_ka/1000
# sigma_esc = calc_xsm_sigma(E_esc)

spec_cu = hypermet_func(E_cu,ene_out_512,a,b,c,i_tail0,i_shelf,alpha,beta,gamma,p,q,r)*total_counts

#Mo
# sigma_main = calc_xsm_sigma(E_mo)
# i_tail = i_tail0*(E_mo/15)**gamma
# i_esc = np.interp(E_mo,esc_peak_int_arr[:,0],esc_peak_int_arr[:,1])
# E_esc = E_mo - si_ka/1000
# sigma_esc = calc_xsm_sigma(E_esc)

spec_mo = hypermet_func(E_mo,ene_out_512,a,b,c,i_tail0,i_shelf,alpha,beta,gamma,p,q,r)*total_counts




####### Escape peak intensity reverse calculations

calc_i_esc = lambda p, energy: p[0]*(np.exp(-p[1]*(energy-1.739))+p[2])

err_esp = lambda p,x,y: calc_i_esc(p,x) - y#p[0]*(np.exp(-p[1]*(x-1.739))+p[2]) - y

p0 = [1.7e-2,0.42,0.04]

p1, success = optimize.leastsq(err_esp, p0[:], args=(esc_peak_int_arr[250:1700,0], esc_peak_int_arr[250:1700,1]))


##### Fitting

p_orig = [a,b,c,i_tail0,i_shelf,alpha,beta,gamma,p,q,r]


def fitfunc(p,ene_out_512): #single ene and spec
    return lambda ene: hypermet_func(ene,ene_out_512,*p)

def errfunc(p,ene_out_512,enes,specs,total_counts): #p - a,b,c,i_tail0,i_shelf,alpha,beta,gamma,p,q,r,  np.array of energies and corresponding monoenergetic normalized spectra
    # errs = np.ravel(hypermet_func(enes,ene_out_512,*p) - specs)
    # errs = np.ravel(fitfunc(p,ene_out_512)(enes) - specs)
    errs = np.zeros(specs.shape)
    for i in range(len(enes)):
        # print(enes[i])
        tmp_err = ((fitfunc(p,ene_out_512)(enes[i]))*total_counts - specs[i,:])/np.sqrt(specs[i,:])
        tmp_err[np.isnan(tmp_err)] = 0
        tmp_err[np.isinf(tmp_err)] = 0
        errs[i,:] = tmp_err
    errs = np.ravel(errs)
    return errs

# def fit_hypermet(enes,specs): # list of energies and corresponding monoenergetic spectra

nIter = 5000
p0 = [3000,300,10,5e-3,1e-4,0.3,0.5,5,2e-2,0.5,0.04]

results = np.zeros((nIter,len(p0)))

for i in range(nIter):
    print(i)
    spec_fe_r = np.random.poisson(spec_fe)
    spec_co_r = np.random.poisson(spec_co)
    spec_cu_r = np.random.poisson(spec_cu)
    spec_mo_r = np.random.poisson(spec_mo)

    # spec_fe_r = spec_fe_r/np.sum(spec_fe_r)
    # spec_co_r = spec_co_r/np.sum(spec_co_r)
    # spec_cu_r = spec_cu_r/np.sum(spec_cu_r)
    # spec_mo_r = spec_mo_r/np.sum(spec_mo_r)


    enes = np.array([E_fe,E_co,E_cu,E_mo])
    specs = np.zeros((len(enes),len(ene_out_512)))

    # specs = np.array([spec_fe_r,spec_co_r,spec_cu_r,spec_mo_r])

    specs[0,:] = spec_fe_r
    specs[1,:] = spec_co_r
    specs[2,:] = spec_cu_r
    specs[3,:] = spec_mo_r



    p1, success = optimize.leastsq(errfunc,p0[:],args=(ene_out_512,enes,specs,total_counts))
    results[i,:] = p1


fig = corner.corner(
    results,
    labels=['a','b','c','I_tail0','I_shelf',r'$\alpha$',r'$\beta$',r'$\gamma$','p','q','r'],
    quantiles=[0.16, 0.5, 0.84],
    show_titles=True,
    title_kwargs={"fontsize": 12},
)

fig.set_size_inches(20,20)

### generating srf files

BCF_DIR = '/home/abhilash/abhilash/SOLEXS/srf/hypermet_srf/hypermet_caldbs/bcf'
SDD_number = 1
rsps_dir = '/media/abhilash/bea6a276-1685-41b3-accd-b394aa7c0539/lalitha/suryahome/suryahome/abhilash/caldbs/hypermet_caldbs/rsps'

ene_in, ene_out, ene_out_512, eff_area, srf, srf_512, rsp, rsp_512 = solexs_caldbgen.calc_caldb(BCF_DIR,SDD_number)

for i in range(nIter):
    srf_hyp = hypermet_func(ene_in.mean(axis=1),ene_out_512,*results[i,:]) 
    rsp = solexs_caldbgen.srf2rsp(srf_hyp, eff_area)
    
    #solexs_caldbgen.write_caldb(ene_in, ene_out, ene_out_512, eff_area, srf_hyp, srf_512, rsp, rsp_512, CPF_DIR, SDD_number)
    rsp_file = rsps_dir+'/solexs_hypermet_resp_'+ str(i).zfill(4)+'.rsp'
    solexs_caldbgen.write_rsp(ene_in, ene_out_512, rsp, rsp_file)

