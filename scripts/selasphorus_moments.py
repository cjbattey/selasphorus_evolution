import moments, sys, os, matplotlib, numpy as np
from moments import Misc,Spectrum,Numerics,Manips,Integration,Demographics1D,Demographics2D
from matplotlib import pyplot as plt
os.chdir("/home/cbattey2/selasphorus/moments/")


fs = moments.Spectrum.from_file("cal_ruf_sas.sfs")
fs=fs.fold()
fs.mask[1,:,:]=True #mask singletons in all populations
fs.mask[:,1,:]=True
fs.mask[:,:,1]=True

ns=fs.sample_sizes

def IM(params, ns):
    nCal,nRuf,nSas,nRS,T_crs,T_rs,m_c_rs,m_rs,m_cr,m_cs = params
    sts = moments.LinearSystem_1D.steady_state_1D(ns[0] + ns[1]+ns[2])
    fs = moments.Spectrum(sts)
    fs = moments.Manips.split_1D_to_2D(fs, ns[0], ns[1]+ns[2])
    fs.integrate([nCal, nRS], T_crs, m = np.array([[0, m_c_rs], [m_c_rs, 0]]))
    fs = moments.Manips.split_2D_to_3D_2(fs,ns[1],ns[2])
    migmat=np.array([[0,m_cr,m_cs],
    				 [m_cr,0,m_rs],
    				 [m_cs,m_rs,0]])
    fs.integrate([nCal,nRuf,nSas],T_rs,m=migmat)
    return  fs
    
upper_bound = [10,10,10,10,10,10,2,2,2,2]
lower_bound = [1e-4,1e-4,1e-4,1e-4,1e-4,1e-4,1e-4,1e-4,1e-4,1e-4]


poptg=[np.random.uniform(lower_bound[x],upper_bound[x]) for x in range(10)]
poptg=moments.Inference.optimize_log(poptg, fs, IM,lower_bound=lower_bound,
                                     upper_bound=upper_bound,verbose=True, 
                                     maxiter=1)

model=IM(poptg, ns)
ll_model=moments.Inference.ll_multinom(model,fs)
theta = moments.Inference.optimal_sfs_scaling(model, fs)
L=950000000
nref=theta/(4*2.3e-9*L)
nCal=poptg[0]*nref
nRuf=poptg[1]*nref
nSas=poptg[2]*nref
nRS=poptg[3]*nref
T_rs=poptg[5]*2*nref
T_rsc=poptg[4]*2*nref+T_rs
m_c_rs=poptg[6]/(2*nref)
m_rs=poptg[7]/(2*nref)
m_cr=poptg[8]/(2*nref)
m_cs=poptg[9]/(2*nref)

out1=[nref,nCal,nRuf,nSas,nRS,T_rs,T_rsc,m_c_rs,m_rs,m_cr,m_cs]
out1="\t".join(map(str,out1))+"\n"
f=open("/home/cbattey2/Dropbox/selasphorus/realparams.txt",'a')
f.write(out1)
f.close()
	
out2="\t".join(map(str,poptg))+"\n"
f=open("/home/cbattey2/Dropbox/selasphorus/modelparams.txt",'a')
f.write(out1)
f.close()

