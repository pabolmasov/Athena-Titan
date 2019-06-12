from numpy import *

import glob

ifplot = True
if(ifplot):
    import matplotlib
    from matplotlib.pyplot import *

def MPItabread(ddir, filemask, fileno, np, ncol = 6):
    # ddir -- directory stem
    # filemask -- file name prefix
    # fileno -- (integer) number of the output file
    # np -- number of the cores involved
    # ncol -- number of the column with the analyzed data
    
    k = arange(np, dtype=integer)
    #    files = ddir + "/id"+str(k)+"/"+filemask+"-id"+str(k)+'{:04d}'.format(fileno)+".tab"

    for kk in arange(np):
        if(kk>0):
            filek = ddir + "/id"+str(kk)+"/"+filemask+"-id"+str(kk)+"."+'{:04d}'.format(fileno)+".tab"
        else:
            filek = ddir + "/id"+str(kk)+"/"+filemask+"."+'{:04d}'.format(fileno)+".tab"
        print(filek)
        lines = loadtxt(filek, comments="#", unpack=False)
        if(kk == 0):
            x = lines[:,2]  ;  y = lines[:,3]  ; q = lines[:, ncol]
        else:
            x = concatenate((x, lines[:,2]))
            y = concatenate((y, lines[:,3]))
            q = concatenate((q, lines[:,ncol]))

    # reshape:
    xun = unique(x) ; yun = unique(y)
    nx = size(xun) ; ny = size(yun)

    # restoring the grid:
    dx = xun[1]-xun[0] ; dy = yun[1]-yun[0] # as the grid is uniform
    kx = rint((x-xun.min())/dx) ;   ky = rint((y-yun.min())/dy)
    i = (kx + ky * nx).astype(int) # effective number of the element
    xnew = zeros(nx*ny) ; ynew = zeros(nx*ny) ; qnew = zeros(nx*ny)
    xnew[i] = x[:] ; ynew[i] = y[:] ; qnew[i] = q[:]
    
    x=reshape(xnew, [nx, ny]) ; y = reshape(ynew, [nx, ny])
    q=reshape(qnew, [nx, ny])

    return x, y, q    
            
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

def xfft(fileno, prefix = "khiso", permode = 4, plots = True, ncol = 6):

    #    x, y, q = tabread("../bin/"+prefix+"{:04d}".format(fileno)+".tab", ncol = 6)
    x, y, q = MPItabread("../bin", prefix, fileno, np = 3, ncol = ncol)
    print(shape(q))
    nx = size(unique(x))
    q_f = fft.rfft(q-q.mean(), axis=0)
    q_sp = (abs(q_f)**2).mean(axis=-1)
    q_dsp = (abs(q_f)**2).std(axis=-1)
    f = fft.rfftfreq(nx, d=(x.max()-x.min())/double(nx))
    
    outfile = "../bin/"+prefix+"{:04d}".format(fileno) # output file prefix
    print(outfile)
    if(plots & ifplot):
        lev1 = q.min() ; lev2 = q.max() ; nlev = 50
        levs = (lev2-lev1) * arange(nlev) / double(nlev-1) + lev1
        clf()
        contourf(x, y, q, levels = levs)
        xlabel('x') ; ylabel('y')
        savefig(outfile+".png")
        clf()
        plot(x[:, 0], q.mean(axis=-1), 'k-')
        plot(x[:, 0], q.std(axis=-1), 'r-')
        savefig(outfile+'_vxcurve.png')
    
        clf()
        plot(f, q_sp, 'k.')
        plot(f*0.+permode, q_sp, 'r-')
        errorbar(f, q_sp, yerr = q_dsp, fmt = 'k.')
        xlabel('$k$') ;   ylabel('PDS')
        xscale('log') # ; yscale('log')
        savefig(outfile+'_vx.png')
    return f, q_sp
        
def fourier_analysis(prefix, nf, permode = 4):
    #    flist = np.sort(glob.glob(prefix+"[0-9][0-9][0-9][0-9].tab"))
        
    for k in arange(nf):
        f, sp = xfft(k, prefix = prefix, permode = permode, plots = False)
        if k == 0:
            nfreq = size(f)
            t2 = zeros([nf, nfreq])
            f2 = zeros([nf, nfreq])
            sp2 = zeros([nf, nfreq])
        f2[k,:] = f ; sp2[k,:] = sp
        t2[k, :] =  double(k) # to be replaced by real time

    if(ifplot):
        clf()
        subplot(211)
        contourf(t2, f2, log(sp2), nlevels=50)
        colorbar()
        yscale('log') ; ylim(1,100)
        xlabel('$t$, s') ; ylabel('$f$, Hz')
        subplot(212)
        plot(t2[:,permode], sp2[:,permode], label="$k = "+str(permode)+"$")
        plot(t2[:,permode*2], sp2[:,permode*2], label="$k = "+str(permode*2)+"$")
        plot(t2[:,1], sp2[:,1], label="$k = "+str(1)+"$")
        plot(t2[:,2], sp2[:,2], label="$k = "+str(2)+"$")
        plot(t2[:,11], sp2[:,11], label="$k = "+str(11)+"$")
        ylim(sp2.max()*1.e-5, sp2.max())
        legend()
        yscale('log') ;  xlabel('$t$, s') ; ylabel('PDS power')
        savefig('fan.png')
        close()
        
