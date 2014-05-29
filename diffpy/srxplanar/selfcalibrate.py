import numpy as np
import scipy as sp
import os
from functools import partial
from scipy.optimize import minimize, leastsq, fmin_bfgs, fmin_l_bfgs_b, fmin_tnc, minimize_scalar, fmin_powell

def halfcut(p, srx, image, xycenter, qind=[50, 1000], show=False, mode='x', output=0):
    if mode == 'x':
        srx.updateConfig(xbeamcenter=p)
    elif mode == 'y':
        srx.updateConfig(ybeamcenter=p)
    elif mode == 'tilt':
        srx.updateConfig(tiltd=p)
    elif mode == 'rotation':
        srx.updateConfig(rotationd=p)
    elif mode == 'all':
        srx.updateConfig(xbeamcenter=p[0],
                         ybeamcenter=p[1],
                         rotationd=p[2],
                         tiltd=p[3])
    kwargs = {'savename':None,
              'savefile':False,
              'flip':False,
              'correction':False,
              }
    if mode != 'y':
        # x half
        mask = np.zeros((srx.config.ydimension, srx.config.xdimension), dtype=bool)
        mask[:, :xycenter[0]] = 1
        res1 = srx.integrate(image, extramask=mask, **kwargs)
        chi1 = res1['chi'][1][qind[0]:qind[1]]
    
        mask = np.logical_not(mask)
        res2 = srx.integrate(image, extramask=mask, **kwargs)
        chi2 = res2['chi'][1][qind[0]:qind[1]]
        
    if mode != 'x':
        # y half
        mask = np.zeros((srx.config.ydimension, srx.config.xdimension), dtype=bool)
        mask[:xycenter[1], :] = 1
        res3 = srx.integrate(image, extramask=mask, **kwargs)
        chi3 = res3['chi'][1][qind[0]:qind[1]]
    
        mask = np.logical_not(mask)
        res4 = srx.integrate(image, extramask=mask, **kwargs)
        chi4 = res4['chi'][1][qind[0]:qind[1]]
        
    if mode == 'x':
        rv = chi1 - chi2
        rv = rv / rv.max()
    elif mode == 'y':
        rv = chi3 - chi4
        rv = rv / rv.max()
    else:
        r1 = chi1 - chi2
        r2 = chi3 - chi4
        # r3 = chi1 - chi3
        # r4 = chi2 - chi4
        # rv = np.concatenate([r1 / r1.max(), r2 / r2.max(), r3 / r3.max(), r4 / r4.max()])
        rv = np.concatenate([r1 / r1.mean(), r2 / r2.mean()])
    
    rv0 = np.sum(rv ** 2)
    print p
    print rv0
    if output == 0:
       rv = rv0

    if show:
        print p
        print rv
        import matplotlib.pyplot as plt
        plt.figure(1)
        plt.clf()
        if mode != 'y':
            plt.plot(res1['chi'][0], res1['chi'][1])
            plt.plot(res2['chi'][0], res2['chi'][1])
        if mode != 'x':
            plt.plot(res3['chi'][0], res3['chi'][1])
            plt.plot(res4['chi'][0], res4['chi'][1])
        plt.show()
    return rv


def selfCalibrateX(srx, image, qmax=20.0, mode='all', output=0):
    bak = {}
    for opt in ['uncertaintyenable', 'integrationspace', 'qmax', 'qstep']:
        bak[opt] = getattr(srx.config, opt)
    
    xycenter = [int(srx.config.xbeamcenter), int(srx.config.ybeamcenter)]
    
    srx.updateConfig(uncertaintyenable=False,
                     integrationspace='qspace',
                     # qmax=qmax,
                     qstep=0.02)
    qind = [50, int(qmax / 0.02)]
    
    srx.prepareCalculation(pic=image)
    srxconfig = srx.config
    image = np.array(srx._getPic(image))
    
    func = partial(halfcut, srx=srx, image=image, qind=qind, mode=mode, output=output,
                   xycenter=xycenter, show=False)

    if mode == 'x':
        p0 = [srxconfig.xbeamcenter]
        bounds = (p0[0] - 3, p0[0] + 3)
    elif mode == 'y':
        p0 = [srxconfig.ybeamcenter]
        bounds = (p0[0] - 3, p0[0] + 3)
    elif mode == 'tilt':
        p0 = [srxconfig.tiltd]
        bounds = (p0[0] - 5, p0[0] + 5)
    elif mode == 'rotation':
        p0 = [srxconfig.rotationd]
        bounds = (0, 360)
    elif mode == 'all':
        p0 = [srxconfig.xbeamcenter, srxconfig.ybeamcenter, srxconfig.rotationd, srxconfig.tiltd]
        bounds = [[p0[0] - 2, p0[0] + 2], [p0[0] - 2, p0[0] + 2], [0, 360], [srxconfig.tiltd - 10, srxconfig.tiltd + 10]]
    
    if output == 0:
        if mode != 'all':
            rv = minimize_scalar(func, bounds=bounds, method='Bounded')
            p = [rv.x]
        else:
            # rv = minimize(func, p0, method='L-BFGS-B', bounds=bounds, options={'xtol':0.001})
            rv = minimize(func, p0, method='Powell', bounds=bounds, options={'xtol':0.001})
            p = rv.x
    else:
        rv = leastsq(func, p0, epsfcn=0.001)
        p = rv[0]
    
    print p
    if mode == 'x':
        srx.updateConfig(xbeamcenter=p[0], **bak)
    elif mode == 'y':
        srx.updateConfig(ybeamcenter=p[0], **bak)
    elif mode == 'tilt':
        srx.updateConfig(tiltd=p[0], ** bak)
    elif mode == 'rotation':
        srx.updateConfig(rotation=p[0], ** bak)
    elif mode == 'all':
        srx.updateConfig(xbeamcenter=p[0], ybeamcenter=p[1], rotationd=p[2], tiltd=p[3], ** bak)
    return p

def selfCalibrate(srx, image, full=True):
    p = []
    if not full:
        p = selfCalibrateX(srx, image, mode='x')
        p = selfCalibrateX(srx, image, mode='y')
        p = selfCalibrateX(srx, image, mode='tilt')
        p = selfCalibrateX(srx, image, mode='rotation')
    else:
        p = selfCalibrateX(srx, image, mode='all')
    return p
