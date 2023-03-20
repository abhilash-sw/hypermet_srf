#####################################################
# @Author: Abhilash Sarwade
# @Date:   2021-11-27 13:56:55
# @email: sarwade@ursc.gov.in
# @File Name: fit_fe55.py
# @Project: hypermet_srf

# @Last Modified time: 2023-03-16 10:28:26
#####################################################

import hypermet_simulation # doesn't work, have to import function only
from scipy import optimize
#from solexs_pipeline.calibration_fit_routines import gaussian, two_gaussian
from solexs_pipeline.calibration_fit_routines import fit_two_gaussian
import numpy as np
from xraydb import XrayDB
import matplotlib.pyplot as plt

xdb = XrayDB()


fe_spec = np.loadtxt('EM_data/long_integration_fe55_spec_solexsEM.dat')
fe_spec = fe_spec.reshape(512,2).sum(axis=1)

ch = np.arange(512)
# calibration from http://cdeg.ursc.dos.gov.in:8081/gitlab/sarwade/SoLEXS_EM_burn_in_20200917/blob/master/ipynb/e_ch_analysis_60eV.ipynb

m1 = 56.73/2 # from EM data
c1 = 16.54 # from EM data

ene = (ch*m1+c1)*1e-3
ene2 = ene + m1*1e-3

ene_out_512 = np.vstack((ene,ene2)).T


### Fitting Gaussians
# def fitfunc_gauss5(p,ch,spec):
p1, p1_err = fit_two_gaussian(ch,fe_spec,guess=[1e7,208,6,2e6,228,12],lower=190,upper=240)

p2, p2_err = fit_two_gaussian(ch,fe_spec,guess=[8e3,158,3,1.1e4,167,3],lower=155,upper=172)

p3, p3_err = fit_two_gaussian(ch,fe_spec,guess=[9e3,130,1,4e4,146,3],lower=124,upper=152)


E_fe_ka = xdb.xray_lines('Mn')['Ka1'].energy/1e3
E_fe_kb = xdb.xray_lines('Mn')['Kb1'].energy/1e3
E_si_ka = xdb.xray_lines('Si')['Ka1'].energy/1e3
E_fe_ka_esc = E_fe_ka - E_si_ka
E_fe_kb_esc = E_fe_kb - E_si_ka

E_ti_ka = xdb.xray_lines('Ti')['Ka1'].energy/1e3
E_unk = 3.7

e_lines = [E_fe_ka, E_fe_kb, E_fe_ka_esc, E_fe_kb_esc, E_ti_ka, E_unk ]

fwhms = np.array([p1[2],p1[5],p3[5],p2[5],p2[2],p3[2]])*m1*2.35

#########
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

E_fe_ka = xdb.xray_lines('Mn')['Ka1'].energy/1e3
E_fe_kb = xdb.xray_lines('Mn')['Kb1'].energy/1e3

ratio_int = xdb.xray_lines('Mn')['Kb1'].intensity/xdb.xray_lines('Mn')['Ka1'].intensity / 0.7

spec_ka = hypermet_func(E_fe_ka,ene_out_512,a,b,c,i_tail0,i_shelf,alpha,beta,gamma,p,q,r)
spec_kb = hypermet_func(E_fe_kb,ene_out_512,a,b,c,i_tail0,i_shelf,alpha,beta,gamma,p,q,r)

spec_ka = spec_ka[0]
spec_kb = spec_kb[0]


def fitfunc(p,ene_out_512,E_fe_ka,E_fe_kb,ratio_int):
    p[0] = 4000
    p[1] = 350
    p[2] = 14
    p[3] = 0.0621
    p[4] = 0.0005

    spec_ka = hypermet_func(E_fe_ka,ene_out_512,*p)
    spec_kb = hypermet_func(E_fe_kb,ene_out_512,*p)

    spec_ka = spec_ka[0]
    spec_kb = spec_kb[0]

    spec_add = spec_ka+ratio_int*spec_kb

    # normal_fe_spec = fe_spec/np.sum(fe_spec)
    normal_spec_add = spec_add/np.sum(spec_add)
    return normal_spec_add

def errfunc(p,fe_spec,ene_out_512,E_fe_ka,E_fe_kb,ratio_int):

    normal_spec_add = fitfunc(p,ene_out_512,E_fe_ka,E_fe_kb,ratio_int)
    err = normal_fe_spec - normal_spec_add

    return err

p_orig = [a,b,c,i_tail0,i_shelf,alpha,beta,gamma,p,q,r]

p0 = [4000,500,20,5e-3,1e-4,0.3,0.5,5,2e-2,0.5,0.04]

p1, success = optimize.leastsq(errfunc,p0[:],args=(fe_spec,ene_out_512,E_fe_ka,E_fe_kb,ratio_int))



normal_fe_spec = fe_spec/np.sum(fe_spec)

plt.style.use(['science','grid','no-latex','notebook'])

fig, ax = plt.subplots()

ax.semilogy(ene_out_512[:,0],normal_fe_spec,label='Fe55 Observation',lw=5)
ax.semilogy(ene_out_512[:,0],fitfunc(p_orig,ene_out_512,E_fe_ka,E_fe_kb,ratio_int),label='Hypermet Fit',lw=3,color='r')
ax.set_ylim([1e-6,1])
ax.set_xlim([0,8])
ax.legend()

ax.set_xlabel('Energy (keV)')
ax.set_ylabel('Normalized Counts')

fig.savefig('fe55_data_hypermet_fit.png')
fig.savefig('fe55_data_hypermet_fit.pdf')





####
ene_peak_fe = [xdb.xray_lines('Mn')['Ka1'].energy,xdb.xray_lines('Mn')['Kb1'].energy,xdb.xray_lines('Mn')['Ka1'].energy-xdb.xray_lines('Si')['Ka1'].energy]


fig, ax = plt.subplots()
# fig.set_size_inches(16,9)

ax.semilogy(ene_out_512[:,0],fe_spec/10,label='Fe55 Observation',lw=5)
ax.semilogy(ene_out_512[:,0],fitfunc(p_orig,ene_out_512,E_fe_ka,E_fe_kb,ratio_int)*np.max(fe_spec),label='Hypermet Fit',lw=3,color='r')
ax.set_ylim([1,1e8])
ax.set_xlim([0,8])
ax.legend(loc=2)

ax.set_xlabel('Energy (keV)')

ax.axvline([0.715],color='k',ls='--',lw=3)
ax.axvline([ene_peak_fe[0]/1000],color='k',ls='--',lw=3,ymax=np.log(3.5e4)/np.log(1e6))
ax.axvline([ene_peak_fe[1]/1000],color='k',ls='--',lw=3,ymax=np.log(1e4)/np.log(1e6))
ax.axvline([ene_peak_fe[2]/1000],color='k',ls='--',lw=3,ymax=np.log(1.5e3/2.5)/np.log(1e6))
ax.axvline([ene_peak_fe[1]/1000 - 1.739],color='k',ls='--',lw=3,ymax=np.log(3.5e2/1.6)/np.log(1e6))

ax.text(ene_peak_fe[0]/1000-0.15,1500000,'5.9 keV\n(FWHM = 192.34 eV)',fontsize=18,ha='right')
ax.text(ene_peak_fe[1]/1000+0.15,400000,r'6.4 keV',fontsize=18,ha='left')
ax.text(ene_peak_fe[2]/1000-0.10,10000,r'4.2 keV',fontsize=18,ha='right')
ax.text(ene_peak_fe[1]/1000+0.00-1.739,5000,r'4.7 keV',fontsize=18,ha='left')
ax.text(0.75,10000,'0.7 keV \n(Low Energy Thresh)',fontsize=18,ha='left')
ax.text(4.99,50000,'Escape Peaks',fontsize=18,ha='right',fontweight='bold')
ax.text(6.8,10000000,r'Fe$^{55}$ X-Rays',fontsize=18,ha='right',fontweight='bold')

fig.savefig('fe55_data_hypermet_fit_an.png')
fig.savefig('fe55_data_hypermet_fit_an.pdf')