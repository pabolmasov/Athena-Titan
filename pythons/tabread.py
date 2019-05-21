from numpy import *
import matplotlib
from matplotlib.pyplot import *

import glob

def tabread(infile, ncol = 6):
    # reading a .tab athena output
    # input: input file and column number for the analysed quantity
    # output: (2D) x, y, and the quantity in question
    lines = loadtxt(infile, comments="#", unpack=False)
    x = lines[:,2]  ;  y = lines[:,3]  ; q = lines[:, ncol]
    # reshape:
    xun = unique(x) ; yun = unique(y)
    nx = size(xun) ; ny = size(yun)

    x=transpose(reshape(x, [ny, nx])) ; y = transpose(reshape(y, [ny, nx]))
    q= transpose(reshape(q, [ny, nx]))

    return x, y, q

def xfft(infile, permode = 4, plots = True):

    x, y, q = tabread(infile, ncol = 6)
    nx = size(unique(x))
    q_f = fft.rfft(q-q.mean(), axis=0)
    q_sp = (abs(q_f)**2).mean(axis=-1)
    q_dsp = (abs(q_f)**2).std(axis=-1)
    f = fft.rfftfreq(nx, d=(x.max()-x.min())/double(nx))

    if(plots):
        clf()
        plot(x[:, 0], q.mean(axis=-1), 'k-')
        plot(x[:, 0], q.std(axis=-1), 'r-')
        savefig(infile+'_vxcurve.png')
    
        clf()
        plot(f, q_sp, 'k.')
        plot(f*0.+permode, q_sp, 'r-')
        errorbar(f, q_sp, yerr = q_dsp, fmt = 'k.')
        xlabel('$k$') ;   ylabel('PDS')
        xscale('log') # ; yscale('log')
        savefig(infile+'_vx.png')
    return f, q_sp
        
def fourier_analysis(prefix, permode = 4, nfilter = 100):
    flist = np.sort(glob.glob(prefix+"[0-9][0-9][0-9][0-9].tab"))
    
    nf = size(flist)
    if(nfilter < nf):
        flist = flist[:nfilter]
        nf = nfilter
    
    for k in arange(nf):
        infile = flist[k]
        print("reading "+str(infile))
        f, sp = xfft(infile, permode = permode, plots = False)
        if k == 0:
            nfreq = size(f)
            t2 = zeros([nf, nfreq])
            f2 = zeros([nf, nfreq])
            sp2 = zeros([nf, nfreq])
        f2[k,:] = f ; sp2[k,:] = sp
        t2[k, :] =  double(k) # to be replaced by real time

    clf()
    subplot(211)
    contourf(t2, f2, log(sp2), nlevels=50)
    colorbar()
    yscale('log') ; ylim(1,100)
    xlabel('$t$, s') ; ylabel('$f$, Hz')
    subplot(212)
    plot(t2[:,permode], sp2[:,permode], label="$k = "+str(permode)+"$")
    plot(t2[:,permode*2], sp2[:,permode*2], label="$k = "+str(permode*2)+"$")
    plot(t2[:,2], sp2[:,2], label="$k = "+str(2)+"$")
    plot(t2[:,11], sp2[:,11], label="$k = "+str(11)+"$")
    ylim(sp2.max()*1.e-5, sp2.max())
    legend()
    yscale('log') ;  xlabel('$t$, s') ; ylabel('PDS power')
    savefig('fan.png')
    close()
