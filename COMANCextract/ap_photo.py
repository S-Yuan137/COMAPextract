import matplotlib.pyplot as plt
from astropy import units as u
from astropy.coordinates import SkyCoord
import numpy as np
import scipy.optimize as opt
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.patches import Ellipse

class apPhotometry:
    '''
    Class to hold photometry tools and utilities.
    '''

    def __init__(self,map,source_coords,method='apPhoto',plotout=False):
        if type(source_coords) is not SkyCoord:
            self.source_coords = SkyCoord(source_coords[0],
                                     source_coords[1],
                                     frame='galactic',
                                     unit=u.degree)
        else:
            self.source_coords = source_coords
        del source_coords

        if method=='mostLikely':
            self.out = self.mostLikely(map,plotout=plotout)
        if method=='apPhoto':
            self.out = self.apPhoto(map,plotout=plotout)

    def apPhoto(self,map,plotout=False):
        print('Aperture Photometry')
        coords = map.coordinates()
        coord_rad = self.source_coords.separation(coords)
        x,y = np.meshgrid(np.arange(map.map.shape[1]),
                          np.arange(map.map.shape[0]))
        x0,y0 = map.wcs.world_to_pixel(self.source_coords)
        coord_theta = np.arctan((x-x0)/(y-y0))
        coord_theta[y-y0<0] += np.pi
        coord_theta[coord_theta<0] += 2*np.pi
        map_med = map.map - np.nanmedian(map.map)

        # FIND AN APPROXIMATE VALUE FOR SIGMA
        r_max,r_step=15,0.1
        flux  = np.zeros(int(r_max/r_step))
        for n in range(int(r_max/r_step)):
            r = np.arange(0,r_max,r_step)[n]
            flux[n] = np.nansum(map_med[coord_rad.arcmin<=r])
        dflux = np.zeros(flux.shape)
        dflux[1:-1] = (flux[2:]-flux[:-2])/(2*r_step)
        X = np.argmax(np.all([dflux<max(dflux)*0.025,
                              flux>max(flux)*0.2],axis=0))
        HWHM = np.argmax(flux>=0.5*flux[X])
        sigma = HWHM / np.sqrt(2*np.log(2))
        sigma_arcmin = sigma*r_step

        # USE THIS TO DEFINE 'ZONES'
        ap_rad = [3,3,3.5]

        coord_rad_sig = (coord_rad.arcmin/sigma_arcmin)
        coord_zone = np.zeros(map.map.shape)
        coord_zone[:] = np.nan
        coord_zone[coord_rad_sig<=ap_rad[0]] = 1
        coord_zone[np.all([ap_rad[1]<coord_rad_sig,
                           coord_rad_sig<=ap_rad[2]],axis=0)] = 2

        in_sum = np.nansum(map.map[coord_zone==1])
        out_sum = np.nansum(map.map[np.all([coord_zone==2],axis=0)])
        in_hit = np.nansum(coord_zone==1)
        out_hit = np.nansum(coord_zone==2)
        S = in_sum - out_sum*(in_hit/out_hit)
        Se = np.nanstd(map.map[coord_zone==2])*np.sqrt(in_hit)

        if plotout==False:
            pass
        else:
            fig = plt.figure(figsize=(8.27,11.69),dpi=70) # A4
            plt.subplots_adjust(hspace=0.27,wspace=0.27)
            plt.suptitle('Most Likely Fit Summary ({} {}GHz)'.format(map.name,map.freq))

            plt.gcf().text(0.1, 0.91,
                           r'$S_\nu =$ {0:.3f} $\pm$ {1:.3f} Jy'.format(S,Se),
                           fontsize=12)

            ax1a = plt.subplot(321,projection=map.wcs)
            ax1a.imshow(map.map,cmap='jet')
            ax1a.set_xlabel(r'$\ell$',labelpad=0.4)
            ax1a.set_ylabel(r'$b$',labelpad=0.6)
            ax1a.set_title('Data')

            ax1b = plt.subplot(322,projection=map.wcs)
            ax1b.imshow(map.map,cmap='Greys')
            ax1b.imshow(coord_zone,cmap='jet',alpha=0.4)
            ax1b.grid(color='black', ls=':')
            ax1b.axvline(x0,color='orange')
            ax1b.axhline(y0,color='green')
            ax1b.set_xlabel(r'$\ell$',labelpad=0.4)
            ax1b.set_ylabel(r'$b$',labelpad=0.6)
            ax1b.set_title('Fit Parameters')

            ax2a = plt.subplot(312)
            ax2a.axvspan(0,sigma_arcmin*ap_rad[0],alpha=0.2,color='blue')
            ax2a.axvspan(sigma_arcmin*ap_rad[1],sigma_arcmin*ap_rad[2],alpha=0.2,color='red')
            ax2a.plot(np.arange(0,r_max,r_step),flux,color='orange',label='Integrated flux')
            ax2a.plot(np.arange(0,r_max,r_step),dflux,color='green',label='Rate of Change')
            for n in range(int(r_max//(4.5/2))+1):
                ax2a.axvline(n*4.5/2,color='black',ls='--',alpha=0.3)
            ax2a.axvline(X*r_step)
            ax2a.axvline(HWHM*r_step)
            ax2a.axvline(sigma*r_step)
            ax2a.set_xlabel(r'$r$ [arcmin]',labelpad=0.4)
            ax2a.set_ylabel(r'$S_\nu$ [Jy/pix]',labelpad=0.6)
            ax2a.set_xlim(0,r_max)
            ax2a.legend()
            ax2a.set_title('Radius Fitting Plot')

            ax3a = plt.subplot(325)
            ax3a.plot((np.arange(map.map.shape[0])-x0)*0.2,
                      map.map[:,int(y0)],c='green')
            ax3a.set_xlabel(r'$x-x_0$ [arcmin]',labelpad=0.4)
            ax3a.set_ylabel(r'$S_\nu$ [Jy/pix]',labelpad=0.6)
            ax3a.set_title('X Axis Cut')

            ax3b = plt.subplot(326)
            ax3b.plot((np.arange(map.map.shape[1])-y0)*0.2,
                      map.map[int(x0),:],c='orange')
            ax3b.set_xlabel(r'$y-y_0$ [arcmin]',labelpad=0.4)
            ax3b.set_ylabel(r'$S_\nu$ [Jy/pix]',labelpad=0.6)
            ax3b.set_title('Y Axis Cut')

            if plotout==True:
                plt.show()
            else:
                plotout.savefig()
                plt.close()
        return S, Se


    def mostLikely(self,map,plotout=False):
        print('Most Likely')
        x0,y0 = map.wcs.world_to_pixel(self.source_coords)
        coords = map.coordinates()
        nll = lambda *args: -self.lnlike(*args)
        p0=[np.nanmax(map.map)-np.nanmedian(map.map), # Amplitude (map_unit)
            self.source_coords.l.degree, # Longitude (degree)
            self.source_coords.b.degree, # Latitude (degree)
            4.5/60, 4.5/60, # S-Maj. and S-Min. sigma (deg)
            np.radians(0), # Angle of rotation (radian)
            np.nanmedian(map.map), # Intercept (map_unit)
            0,0, # Slopes in x and y (map_unit/deg)
            0.0005] # Estimated noise (map_unit)
        eps = 1e-10
        result = opt.minimize(nll,p0,
                              args=(coords.l.degree,
                                    coords.b.degree,
                                    map.map),
                              method='L-BFGS-B',
                              options={'ftol':1.e-16}
                              )

        ftol = 1.e-16
        tmp_i = np.zeros(len(result.x))
        uncert = np.zeros(len(result.x))
        for i in range(len(result.x)):
            tmp_i[i] = 1.0
            hess_inv_i = result.hess_inv(tmp_i)[i]
            uncertainty_i = np.sqrt(max(1, abs(result.fun)) * ftol * hess_inv_i)
            tmp_i[i] = 0.0
            uncert[i] = uncertainty_i
            # print('x^{0} = {1:12.4e} ± {2:.1e}'.format(i, result.x[i], uncertainty_i))

        S = 2*np.pi*result['x'][0]*result['x'][3]*result['x'][4]*3600/0.04
        Se = S*np.sqrt( (uncert[0]/result['x'][0])**2 +
                          (uncert[3]/result['x'][3])**2 +
                          (uncert[4]/result['x'][4])**2  )
        if plotout==False:
            pass
        else:
            fig = plt.figure(figsize=(8.27,11.69),dpi=70) # A4
            plt.subplots_adjust(hspace=0.27,wspace=0.27)
            plt.suptitle('Aperture Photometry Summary ({} {}GHz)'.format(map.name,map.freq))

            plt.gcf().text(0.1, 0.91,
                           r'$S_\nu =$ {0:.3f} $\pm$ {1:.3f} Jy'.format(S,Se),
                           fontsize=12)

            Min,Max = np.nanmin(map.map),np.nanmax(map.map)

            ax1a = plt.subplot(321,projection=map.wcs)
            ax1a.imshow(map.map,cmap='jet',vmin=Min,vmax=Max)
            ax1a.set_xlabel(r'$\ell$',labelpad=0.4)
            ax1a.set_ylabel(r'$b$',labelpad=0.6)
            ax1a.set_title('Data')

            ax1b = plt.subplot(322,projection=map.wcs)
            ax1b.imshow(map.map,cmap='Greys',vmin=Min,vmax=Max)
            # ax1b.imshow(coord_zone,cmap='jet',alpha=0.4)
            ax1b.grid(color='black', ls=':')
            ax1b.axvline(x0,color='orange')
            ax1b.axhline(y0,color='green')

            print(result['x'])

            Ell = Ellipse((result['x'][1],result['x'][2]),
                          result['x'][3]*2,
                          result['x'][4]*2,
                          result['x'][5],
                          edgecolor='red',
                          facecolor='none',
                          transform=ax1b.get_transform('world'))
            ax1b.add_patch(Ell)

            ax1b.set_xlabel(r'$\ell$',labelpad=0.4)
            ax1b.set_ylabel(r'$b$',labelpad=0.6)
            ax1b.set_title('Fit Parameters')

            ax2a = plt.subplot(323,projection=map.wcs)
            ax2a.imshow(Gauss2D(result['x'],coords.l.degree,coords.b.degree),
                        cmap='jet',vmin=Min,vmax=Max)
            ax2a.set_xlabel(r'$\ell$',labelpad=0.4)
            ax2a.set_ylabel(r'$b$',labelpad=0.6)
            ax2a.set_title('Model')

            ax2b = plt.subplot(324,projection=map.wcs)
            ax2b.imshow(map.map-Gauss2D(result['x'],coords.l.degree,coords.b.degree),
                        cmap='jet',vmin=Min,vmax=Max)
            ax2b.set_xlabel(r'$\ell$',labelpad=0.4)
            ax2b.set_ylabel(r'$b$',labelpad=0.6)
            ax2b.set_title('Data-Model')

            ax3a = plt.subplot(325)
            ax3a.plot((np.arange(map.map.shape[0])-x0)*0.2,
                      Gauss2D(result['x'],
                              coords.l.degree,
                              coords.b.degree)[:,int(y0)],
                      'k-',alpha=0.5)
            ax3a.plot((np.arange(map.map.shape[0])-x0)*0.2,
                      Gauss2D(result['x'],
                              coords.l.degree,
                              coords.b.degree)[:,int(y0)] + result['x'][-1],
                      'k--',alpha=0.5)
            ax3a.plot((np.arange(map.map.shape[0])-x0)*0.2,
                      Gauss2D(result['x'],
                              coords.l.degree,
                              coords.b.degree)[:,int(y0)] - result['x'][-1],
                      'k--',alpha=0.5)
            ax3a.plot((np.arange(map.map.shape[0])-x0)*0.2,
                      map.map[:,int(y0)],c='green')
            ax3a.set_xlabel(r'$x-x_0$ [arcmin]',labelpad=0.4)
            ax3a.set_ylabel(r'$S_\nu$ [Jy/pix]',labelpad=0.6)
            ax3a.set_title('X Axis Cut')

            ax3b = plt.subplot(326)
            ax3b.plot((np.arange(map.map.shape[0])-x0)*0.2,
                      Gauss2D(result['x'],
                              coords.l.degree,
                              coords.b.degree)[int(x0),:],
                      'k-',alpha=0.5)
            ax3b.plot((np.arange(map.map.shape[0])-x0)*0.2,
                      Gauss2D(result['x'],
                              coords.l.degree,
                              coords.b.degree)[int(x0),:] + result['x'][-1],
                      'k--',alpha=0.5)
            ax3b.plot((np.arange(map.map.shape[0])-x0)*0.2,
                      Gauss2D(result['x'],
                              coords.l.degree,
                              coords.b.degree)[int(x0),:] - result['x'][-1],
                      'k--',alpha=0.5)
            ax3b.plot((np.arange(map.map.shape[1])-y0)*0.2,
                      map.map[int(x0),:],c='orange')
            ax3b.set_xlabel(r'$y-y_0$ [arcmin]',labelpad=0.4)
            ax3b.set_ylabel(r'$S_\nu$ [Jy/pix]',labelpad=0.6)
            ax3b.set_title('Y Axis Cut')

            if plotout==True:
                plt.show()
            else:
                plotout.savefig()
                plt.close()
        return S, Se

    def lnlike(self,p,x,y,z):
        amplitude, x_mean, y_mean, x_stddev, y_stddev, intercept, mx, my, theta, eps = p
        model = Gauss2D(p,x,y)
        denom = np.power(eps,2)
        lp = -0.5*np.nansum(np.power((z-model),2)/denom + np.log(denom) + np.log(2*np.pi))
        return lp

def Gauss2D(p,x,y):
    amplitude, x_mean, y_mean, x_stddev, y_stddev, theta, intercept, mx, my, eps = p
    cost2 = np.cos(theta) ** 2
    sint2 = np.sin(theta) ** 2
    sin2t = np.sin(2. * theta)
    xstd2 = x_stddev ** 2
    ystd2 = y_stddev ** 2
    xdiff = x - x_mean
    ydiff = y - y_mean
    a = 0.5 * ((cost2 / xstd2) + (sint2 / ystd2))
    b = 0.5 * ((sin2t / xstd2) - (sin2t / ystd2))
    c = 0.5 * ((sint2 / xstd2) + (cost2 / ystd2))
    bg = mx*x + my*y + intercept
    return amplitude * np.exp(-((a * xdiff ** 2) + (b * xdiff * ydiff) +
                                (c * ydiff ** 2))) + bg
